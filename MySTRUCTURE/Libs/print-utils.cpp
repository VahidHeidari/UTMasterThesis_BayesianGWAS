#include "print-utils.h"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "allele-frequencies.h"
#include "logger.h"
#include "no-admix-params.h"

using namespace std::chrono;

static high_resolution_clock::time_point start_time;

void PrintIterationsInfo(int iter, const Params& params)
{
	if (iter == 0)
		start_time = std::chrono::high_resolution_clock::now();

	const int ITERATIONS = params.GetNumIterations() + params.GetNumBurnins();
	if (iter != 0 && (iter + 1) % UPDATE_FREQ != 0 && iter < ITERATIONS)
		return;

	std::cout << "Iteration #" << iter + 1 << "  of " << ITERATIONS
			<< " (" << (iter < params.GetNumBurnins() ? "BURNIN" : "SAMPLING") << ')';

	if (iter != 0) {
		// Calculate left time.
		const auto end_time = high_resolution_clock::now();
		const double one_loop_dur = duration<double>(end_time - start_time).count() / (double)UPDATE_FREQ;
		//std::cout << "   One loop:"<< one_loop_dur << "    ";
		const double rem_iters = ITERATIONS - iter;
		const double dur = one_loop_dur * rem_iters;					// Total duration
		const int sec = (int)dur % 60;
		const int min = (int)(dur / 60.0) % 60;
		const int hur = (int)(dur / 60.0 / 60.0) % 60;
		std::cout << "    left time = " << hur << ':' << min << ':' << sec << "        ";
		start_time = end_time;
	}

	if (iter == ITERATIONS)
		std::cout << std::endl;
	else
		std::cout << '\r';
}

void PrintResults(const Params& params, const IndivClusters& Z, const AlleleFrequencies& P)
{
	std::ofstream out_props("no_admix_props.txt");
	out_props << std::setprecision(2) << std::fixed;
	out_props << "Clusters:" << std::endl;
	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		out_props << '#' << (i + 1) << "  " << Z[i] << "     ";
		for (int k = 0; k < params.GetNumClusters(); ++k) {
			const double LOG_PROB = LogIndivProb(i, k, params, P);
			out_props << LOG_PROB;
			if (k + 1 < params.GetNumClusters())
				out_props << "  ";
		}
		out_props << std::endl;
	}

	out_props << std::endl;
	out_props << "Allele frequencies:" << std::endl;
	for (int k = 0; k < params.GetNumClusters(); ++k) {
		out_props << '#' << k + 1 << std::endl;
		for (int l = 0; l < params.GetNumLoci(); ++l) {
			out_props << "    ";
			for (int a = 0; a < params.GetNumAlleles(l); ++a)
				out_props << P.GetAlleleFreq(k, l, a) << "    ";
			out_props << std::endl;
		}
		out_props << std::endl;
	}
}

void PrintAdmixResults(const Params& params, const AlleleFrequencies& P, const AdmixProportions& Q)
{
	const auto prec = std::cout.precision();
	logger << std::setprecision(2);
	logger << std::endl;
	logger << "Cluster proportions (Q):" << std::endl;
	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		logger << '#' << i + 1 << "    ";
		for (int k = 0; k < params.GetNumClusters(); ++k) {
			logger << Q.GetOriginProportion(i, k);
			if (k + 1 < params.GetNumClusters())
				logger << "    ";
		}
		logger << std::endl;
	}

	logger << std::endl;
	logger << "Allele frequencies (P):" << std::endl;
	for (int k = 0; k < params.GetNumClusters(); ++k) {
		logger << '#' << k + 1 << std::endl;
		for (int l = 0; l < params.GetNumLoci(); ++l) {
			logger << "    ";
			for (int a = 0; a < params.GetNumAlleles(l); ++a) {
				logger << P.GetAlleleFreq(k, l, a);
				if (a + 1 < params.GetNumAlleles(l))
					logger << "  ";
			}
			logger << std::endl;
		}
	}
	logger << std::setprecision(prec);
}

static void LogFullProbs(const DiseaseSite& site, const DiseaseModel& M, int locus)
{
	logger << "      LOCUS                     : " << locus << std::endl;
	logger << "        START LOCI              : " << site.loci_start << std::endl;

	logger << "        I-TH ALLELE COMBINATION : " << site.i_th_allele_comb << "    ";
	const auto& allele_comb = M.GetAlleleCombination(site.i_th_allele_comb);
	for (int i : allele_comb)
		logger << i << "  ";
	logger<< std::endl;

	logger << "        I-TH LOCI COMBINATION   : " << site.i_th_loci_comb << "    ";
	const auto& loci_comb = M.GetSiteCombination(site.i_th_loci_comb);
	for (int i : loci_comb)
		logger << i + 1 + site.loci_start << "  ";
	logger << std::endl;

	logger << "        OMEGA                   : " << site.omega << std::endl;
	logger << "        W                       : " << site.w << std::endl;
	logger << "        #CASES                  : " << site.num_cases << std::endl;
	logger << "        #CONTROLS               : " << site.num_contorls << std::endl;
	logger << "        PROB                    : " << site.prob << std::endl;
	logger << std::endl;
}

static void LogSimpleProbs(const DiseaseSite& site, int locus)
{
	const auto prec = std::cout.precision();
	logger << std::setprecision(2);
	logger << std::right << std::setw(3) << locus + 1 << ':';
	logger << std::left << std::setw(5) << site.prob << ' ';
	logger << std::setprecision(prec);
}

static void LogSingleProbs(const Params& params, const DiseaseModel& M, bool is_log_full)
{
	for (int k = 0; k < params.GetNumClusters(); ++k) {
		logger << "    #" << k + 1 << std::endl;
		for (int l = 0; l < params.GetNumLoci(); ++l) {
			const DiseaseSite& site = M.GetDiseaseSite(k, l);

			if (is_log_full)
				LogFullProbs(site, M, l);
			else
				LogSimpleProbs(site, l);
		}
		logger << std::endl;
	}
}

void PrintExhMotahariResults(const Params& params, const IndivClusters& Z, const AlleleFrequencies& P,
		const DiseaseModel& M, bool is_log_full)
{
	PrintResults(params, Z, P);

	logger << "Disease Model:" << std::endl;
	if (LOG_BOTH_SITE_PROBS) {
		LogSingleProbs(params, M, true);
		LogSingleProbs(params, M, false);
	} else
		LogSingleProbs(params, M, is_log_full);
}

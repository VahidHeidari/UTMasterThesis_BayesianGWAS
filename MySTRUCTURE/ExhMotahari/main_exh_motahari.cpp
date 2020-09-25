#include <cmath>
#include <vector>

#include "logger.h"
#include "no-admix-params.h"
#include "params.h"
#include "print-utils.h"

static const int epistasis_radius = 10;			/// Radius of epistatic neighborhood loci
static const int avg_diseases = 2;				/// Average causal disease sites in epistasis radius

static Params params;

static IndivClusters Z;							/// Membrance of each individual to clusters
static AlleleFrequencies P;						/// Allele frequencies
static AllelesCounts allele_count;				/// Number of each alleles at each locus
static DiseaseModel M;							/// Disease model for each cluster

static void InitializeExhMotahariParameters()
{
	InitializeNoAdmixParameters(params, Z, P, allele_count);

#ifdef TEST_ENTROPY
	for (unsigned i = 0; i < Z.size(); ++i)
		Z[i] = i < Z.size() / 2 ? 0 : 1;
#endif

	M.Init(params.GetNumIndividuals(), params.GetNumLoci(), params.GetNumClusters(),
		params.GetNumUniqueAlleles(), avg_diseases, epistasis_radius);
}

static bool IsCombEqual(int start_locus, const Chromosome& chrom,
		const std::vector<int>& alleles_comb, const std::vector<int>& loci)
{
	for (unsigned l = 0; l < loci.size(); ++l) {
		const int locus = start_locus + loci[l];
		const int indiv_allele = chrom[locus];
		const int comb_allele = alleles_comb[l];
		if (indiv_allele != comb_allele)
			return false;
	}

	return true;
}

static void CountCombinations2(int start_locus, const std::vector<int>& alleles_comb,
	const std::vector<int>& loci, int k, int& omega, int& w, int& num_cases, int& num_controls)
{
	omega = w = num_cases = num_controls = 0;

	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		if (Z[i] != k)
			continue;

		// Update case and controls counters.
		if (params.IsCaseIndividual(i))
			++num_cases;
		else
			++num_controls;

		// Update infected and uninfected combinations counters.
		for (int c = 0; c < params.GetNumChromosomes(); ++c) {
			const Chromosome& chrom = params.GetChromosome(i, c);
			if (IsCombEqual(start_locus, chrom, alleles_comb, loci)) {
				if (params.IsCaseIndividual(i))
					++w;
				++omega;
			}
		}
	}
}

static void CountCombinations(int start_loci, int i_th_allele_comb,
	int i_th_loci_comb, int k, int& omega, int& w, int& num_cases, int& num_controls)
{
	const AllelesComb& allele_comb = M.GetAlleleCombination(i_th_allele_comb);
	const Permuts& loci_comb = M.GetSiteCombination(i_th_loci_comb);
	CountCombinations2(start_loci, allele_comb, loci_comb, k, omega, w, num_cases, num_controls);
}

static bool IsNaN(double d)
{
	return d != d;
}

static void UpdateM()
{
	M.ClearSites();

	const int LAST_LOCUS = params.GetNumLoci() - M.GetEpistasisRadius() + 1;
	for (int l = 0; l < LAST_LOCUS; ++l) {
		for (int a_comb = 0; a_comb < M.GetNumAllelesCombinations(); ++a_comb) {
			for (int s_comb = 0; s_comb < M.GetNumSitesCombinations(); ++s_comb) {
				for (int k = 0; k < params.GetNumClusters(); ++k) {
					int omega, w, num_cases, num_controls;
					CountCombinations(l, a_comb, s_comb, k, omega, w, num_cases, num_controls);
					M.Update(k, l, a_comb, s_comb, omega, w, num_cases, num_controls);
#ifdef TEST_ENTROPY
					if (omega == 0)
						continue;
					const double INF_RATE = (double)w / (double)omega;
					const double PR = std::log(INF_RATE);
					const double PR_COMP = std::log(1.0 - PR);
					const double Pr = std::exp(num_cases * PR + num_controls * PR_COMP);
					if (IsNaN(Pr) || Pr == 0 || Pr == 1 || Pr < 0.7)
						continue;
					logger << "L:" << l << "    O:" << omega << "    W:" << w << "   Pr:" << Pr;
					logger << "   A:" << a_comb << "  ";
					for (const auto& a : M.GetAlleleCombination(a_comb))
						logger << a + 1 << ' ';
					logger << "    S:" << s_comb << "  ";
					for (const auto& s : M.GetSiteCombination(s_comb))
						logger << s + l + 1 << ' ';
					logger << std::endl;
#endif
				}		// Clusters
			}		// Site combinations
		}		// Allele combinations
	}		// Loci
#ifdef TEST_ENTROPY
		exit(1);
#endif
}

int main()
{
	logger << std::endl << std::endl;							// Log start of application.
	logger << "----- EXHAUSTIVE MOTAHARI -----" << std::endl;
	logger << "Start : " << Time << std::endl;

	if (!params.InitExhMotahari(DEF_EXH_MOTAHARI))				// Read input config.
		return 1;

	params.Print();												// Print configuration parameters.
	InitializeExhMotahariParameters();							// Initialize parameters randomly.

	logger << "Start of MCMC loop:" << std::endl;
	const int ITERATIONS = params.GetNumIterations() + params.GetNumBurnins();
	for (int i = 0; i < ITERATIONS; ++i) {
		UpdateP(params, Z, allele_count, P);
		UpdateZ(params, Z, P);
		UpdateM();
		PrintIterationsInfo(i, params);
	}

	PrintIterationsInfo(ITERATIONS, params);
	PrintExhMotahariResults(params, Z, P, M);					// Print and dump results.
	logger << "End   : " << Time << std::endl;
	return 0;
}

#include "no-admix-params.h"

#include <random>

#include "params.h"

void InitializeNoAdmixParameters(const Params& params, IndivClusters& Z, AlleleFrequencies& P, AllelesCounts& allele_count)
{
	Z.clear();
	Z.resize(params.GetNumIndividuals());
	std::random_device dev;
	std::uniform_int<int> unif(0, params.GetNumClusters() - 1);
	for (int i = 0; i < params.GetNumIndividuals(); ++i)
		Z[i] = unif(dev);

	allele_count.Init(params.GetNumClusters(), params.GetNumLoci(), params.GetNumAlleles());
	P.Init(params.GetNumClusters(), params.GetNumLoci(), params.GetNumAlleles());
}

void CountAlleles(const Params& params, const IndivClusters& Z, AllelesCounts& allele_count)
{
	allele_count.ZeroAll();

	// Count alleles.
	for (int l = 0; l < params.GetNumLoci(); ++l) {
		for (int i = 0; i < params.GetNumIndividuals(); ++i) {
			const int cluster = Z[i];
			for (int c = 0; c < params.GetNumChromosomes(); ++c) {
				const int allele = params.GetAllele(i, c, l);
				allele_count.Increment(cluster, l, allele);
			}
		}
	}
}

void SimulateP(const Params& params, const AllelesCounts& allele_count, AlleleFrequencies& P)
{
	for (int l = 0; l < params.GetNumLoci(); ++l) {
		for (int cluster = 0; cluster < params.GetNumClusters(); ++cluster) {
			const int NUM_PARAMETERS = params.GetNumAlleles(l);
			std::vector<double> parameters(NUM_PARAMETERS);
			for (int a = 0; a < NUM_PARAMETERS; ++a) {
				const int COUNT = allele_count.GetAlleleCount(cluster, l, a);
				parameters[a] = LAMBDA + COUNT;
			}
			SimulateDirichlet(parameters, P.GetAlleleFreqs(cluster, l));
		}
	}
}

void UpdateP(const Params& params, const IndivClusters& Z, AllelesCounts& allele_count, AlleleFrequencies& P)
{
	CountAlleles(params, Z, allele_count);
	SimulateP(params, allele_count, P);
}

double LogIndivProb(int i, int k, const Params& params, const AlleleFrequencies& P)
{
	double log_pr = 0.0;
	for (int l = 0; l < params.GetNumLoci(); ++l) {
		for (int c = 0; c < params.GetNumChromosomes(); ++c) {
			const int allele = params.GetAllele(i, c, l);
			const double p = P.GetAlleleFreq(k, l, allele);
			log_pr += std::log(p);
		}
	}
	return log_pr;
}

static double LogAllIndivProb(int i, std::vector<double>& probs, const Params& params, const AlleleFrequencies& P)
{
	double sum = 0;
	for (int k = 0; k < params.GetNumClusters(); ++k) {
		const double pr = LogIndivProb(i, k, params, P);
		probs[k] = std::exp(pr);
		sum += probs[k];
	}
	const double res = std::log(sum);
	return res;
}

void UpdateZ(const Params& params, IndivClusters& Z, const AlleleFrequencies& P)
{
	std::vector<double> probs(params.GetNumClusters());

	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		for (int k = 0; k < params.GetNumClusters(); ++k) {
			const double indiv_prob = LogIndivProb(i, k, params, P);
			probs[k] = indiv_prob;
		}
		const int NEW_K = GetMaxProbIndex(probs);
		Z[i] = NEW_K;
	}
}

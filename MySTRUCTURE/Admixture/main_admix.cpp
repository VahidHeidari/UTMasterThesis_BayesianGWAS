#include <vector>

#include "allele-frequencies.h"
#include "dists.h"
#include "logger.h"
#include "params.h"
#include "print-utils.h"

static Params params;

static AdmixZ Z;						/// Membrance of each alleles of individual to clusters
static AlleleFrequencies P;				/// Allele frequencies
static AdmixProportions Q;				/// Admixture proprotions of each individual
static AllelesCounts allele_count;		/// Number of each alleles at each locus

static void InitializeParameters()
{
	Z.Init(params.GetNumIndividuals(), params.GetNumChromosomes(), params.GetNumLoci(), params.GetNumClusters());
	P.Init(params.GetNumClusters(), params.GetNumLoci(), params.GetNumAlleles());
	Q.Init(params.GetNumIndividuals(), params.GetNumClusters());
	allele_count.Init(params.GetNumClusters(), params.GetNumLoci(), params.GetNumAlleles());
}

static void CountAllels()
{
	allele_count.ZeroAll();
	for (int l = 0; l < params.GetNumLoci(); ++l) {
		for (int i = 0; i < params.GetNumIndividuals(); ++i) {
			for (int c = 0; c < params.GetNumChromosomes(); ++c) {
				const int cluster = Z.GetOrigin(i, c, l);
				const int allele = params.GetAllele(i, c, l);
				allele_count.Increment(cluster, l, allele);
			}
		}
	}
}

static void SimulateP()
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

static void UpdateP()
{
	CountAllels();
	SimulateP();
}

static void GetZCount(int i, std::vector<int>& origin_count)
{
	for (unsigned k = 0; k < origin_count.size(); ++k)
		origin_count[k] = 0;

	for (int k = 0; k < params.GetNumClusters(); ++k) {
		for (int c = 0; c < params.GetNumChromosomes(); ++c) {
			for (int l = 0; l < params.GetNumLoci(); ++l) {
				const int ORIGIN = Z.GetOrigin(i, c, l);
				++origin_count[ORIGIN];
			}
		}
	}
}

static void UpdateQ()
{
	std::vector<int> origin_count(params.GetNumClusters());
	std::vector<double> parameters(params.GetNumClusters());

	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		GetZCount(i, origin_count);
		for (int k = 0; k < params.GetNumClusters(); ++k)
			parameters[k] = ALPHA + origin_count[k];
		SimulateDirichlet(parameters, Q.GetOriginProportions(i));
	}
}

static void UpdateZ()
{
	std::vector<double> probs(params.GetNumClusters());

	for (int i = 0; i < params.GetNumIndividuals(); ++i) {
		for (int c = 0; c < params.GetNumChromosomes(); ++c) {
			for (int l = 0; l < params.GetNumLoci(); ++l) {
				const int allele = params.GetAllele(i, c, l);
				double sum_probs = 0.0;
				for (int k = 0; k < params.GetNumClusters(); ++k) {
					const double q = Q.GetOriginProportion(i, k);
					const double p = P.GetAlleleFreq(k, l, allele);
					probs[k] = q * p;
					sum_probs += probs[k];
				}
				const int NEW_K = SimulateRouletteWheel(probs, sum_probs);
				Z.SetOrigin(i, c, l, NEW_K);
			}
		}
	}
}

static void UpdateAlpha()
{
	// TODO: Not implemented yet! It is assumed fixed and uniform.
}

int main()
{
	logger << std::endl << std::endl;							// Log start of application.
	logger << "***** ADMIXTURE *****" << std::endl;
	logger << "Start : " << Time << std::endl;

	if (!params.Init(DEF_ADMIX_CONF))							// Read input config.
		return 1;

	params.Print();												// Print configuration parameters.
	InitializeParameters();										// Initialize parameters randomly.

	logger << "Start of MCMC loop:" << std::endl;
	const int ITERATIONS = params.GetNumIterations() + params.GetNumBurnins();
	for (int i = 0; i < ITERATIONS; ++i) {
		UpdateP();
		UpdateQ();
		UpdateZ();
		UpdateAlpha();
		PrintIterationsInfo(i, params);
	}

	PrintIterationsInfo(ITERATIONS, params);
	PrintAdmixResults(params, P, Q);							// Print and dump results.
	logger << "End   : " << Time << std::endl;
	return 0;
}

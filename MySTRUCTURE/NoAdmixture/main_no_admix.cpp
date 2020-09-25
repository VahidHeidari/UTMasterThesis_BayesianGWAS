#include <vector>

#include "allele-frequencies.h"
#include "logger.h"
#include "no-admix-params.h"
#include "params.h"
#include "print-utils.h"

static Params params;

static IndivClusters Z;					/// Membrance of each individual to clusters
static AlleleFrequencies P;				/// Allele frequencies
static AllelesCounts allele_count;		/// Number of each alleles at each locus

int main(int argc, char** argv)
{
	logger << std::endl << std::endl;							// Log start of application.
	logger << "----- NO ADMIXTURE -----" << std::endl;
	logger << "Start : " << Time << std::endl;

	const char* INPUT_FILE = argc < 2 ? DEF_NO_ADMIX_CONF : argv[1];
	logger << "Input config file: " << INPUT_FILE << std::endl;
	if (!params.Init(INPUT_FILE))								// Read input config.
		return 1;

	params.Print();												// Print configuration parameters.
	InitializeNoAdmixParameters(params, Z, P, allele_count);	// Initialize parameters randomly.

	logger << "Start of MCMC loop:" << std::endl;
	const int ITERATIONS = params.GetNumIterations() + params.GetNumBurnins();
	for (int i = 0; i < ITERATIONS; ++i) {
		UpdateP(params, Z, allele_count, P);
		UpdateZ(params, Z, P);
		PrintIterationsInfo(i, params);
	}

	PrintIterationsInfo(ITERATIONS, params);
	PrintResults(params, Z, P);									// Print and dump results.
	logger << "End   : " << Time << std::endl;
	return 0;
}

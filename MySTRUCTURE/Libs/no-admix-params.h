#ifndef NO_ADMIX_PARAMS_H_
#define NO_ADMIX_PARAMS_H_

#include <vector>

#include "allele-frequencies.h"

class Params;

typedef std::vector<int> IndivClusters;		/// Membrance of each individual to clusters

void InitializeNoAdmixParameters(const Params& params, IndivClusters& Z, AlleleFrequencies& P, AllelesCounts& allele_count);
void CountAlleles(const Params& params, const IndivClusters& Z, AllelesCounts& allele_count);
void SimulateP(const Params& params, const AllelesCounts& allele_count, AlleleFrequencies& P);
void UpdateP(const Params& params, const IndivClusters& Z, AllelesCounts& allele_count, AlleleFrequencies& P);
double LogIndivProb(int i, int k, const Params& params, const AlleleFrequencies& P);
void UpdateZ(const Params& params, IndivClusters& Z, const AlleleFrequencies& P);

#endif

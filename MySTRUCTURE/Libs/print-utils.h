#ifndef PRINT_UTILS_H_
#define PRINT_UTILS_H_

#include <vector>

#include "no-admix-params.h"
#include "params.h"

class AdmixProportions;
class AlleleFrequencies;
class DiseaseModel;
class Params;

void PrintIterationsInfo(int iter, const Params& params);
void PrintResults(const Params& params, const IndivClusters& Z, const AlleleFrequencies& P);
void PrintAdmixResults(const Params& params, const AlleleFrequencies& P, const AdmixProportions& Q);
void PrintExhMotahariResults(const Params& params, const IndivClusters& Z, const AlleleFrequencies& P,
	const DiseaseModel& M, bool is_log_full = LOG_FULL_SITE_PROBS);

#endif

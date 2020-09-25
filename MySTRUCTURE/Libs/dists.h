#ifndef DISTS_H_
#define DISTS_H_

#include <vector>

typedef std::vector<int> Permuts;
typedef std::vector<Permuts> PermutsVect;

void SimulateDirichlet(const std::vector<double>& alphas, std::vector<double>& output);
int SimulateRouletteWheel(const std::vector<double>& probs, double sum_probs);
int GetMaxProbIndex(const std::vector<double>& probs);

double PoissonDensity(double x, double mu);

void MakePermutations(int n, PermutsVect& output);
void MakeCombinations(int n, int r, PermutsVect& output);
void FastMakeCombinations(int n, int r, PermutsVect& output_combs);

bool DumpPermutations(const char* file_name, int n, const PermutsVect& perms);
bool DumpCombinations(const char* file_name, int n, int r, const PermutsVect& combs);

#endif

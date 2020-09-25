#include "dists.h"

#include <algorithm>
#include <fstream>
#include <random>

#include "logger.h"

void SimulateDirichlet(const std::vector<double>& alphas, std::vector<double>& output)
{
	std::random_device dev;
	double sum = 0.0;
	const int ALPHAS_SIZE = static_cast<int>(alphas.size());
	for (int i = 0; i < ALPHAS_SIZE; ++i) {
		const double alpha = alphas[i];
		std::gamma_distribution<double> gamma(alpha, 1);
		output[i] = gamma(dev);
		sum += output[i];
	}
	for (int i = 0; i < ALPHAS_SIZE; ++i)
		output[i] /= sum;
}

int SimulateRouletteWheel(const std::vector<double>& probs, double sum_probs)
{
	// Roulette wheel random selection.
	std::random_device dev;
	std::uniform_real<double> unif(0, sum_probs);
	const double u = unif(dev);
	double sum = 0;
	for (int k = 0; k < static_cast<int>(probs.size()); ++k) {
		sum += probs[k];
		if (sum < u)
			continue;

		return static_cast<int>(k);
	}

	return static_cast<int>(probs.size()) - 1;
}

int GetMaxProbIndex(const std::vector<double>& probs)
{
	double max = -std::numeric_limits<double>::max();
	int idx = 0;
	for (int k = 0; k < static_cast<int>(probs.size()); ++k) {
		if (probs[k] > max) {
			max = probs[k];
			idx = k;
		}
	}
	return idx;
}

static double Factorial(int n)
{
	if (n == 0 || n == 1)
		return 1.0;

	if (n < 0) {
		logger << warning << "Factorial(" << n << ") should be greater than zero!" << std::endl;
		return 0;
	}

	double res = 1.0;
	for (int i = 2; i < n; ++i)
		res *= i;
	return res;
}

double PoissonDensity(int k, double lambda)
{
	double res = std::exp(-lambda) * std::pow(lambda, k) / Factorial(k);
	return res;
}

void MakePermutations(int n, PermutsVect& output)
{
	output.clear();

	Permuts perm(n);
	for (int i = 0; i < n; ++i)
		perm[i] = i;

	do {
		output.push_back(perm);
	} while (std::next_permutation(perm.begin(), perm.end()));
}

static bool IsCombDuplicated(int r, const Permuts& perm, const PermutsVect& output)
{
	if (!output.size())
		return false;

	Permuts comb = perm;
	comb.resize(r);
	std::sort(comb.begin(), comb.end());
	for (const auto& o_comb : output)
		if (comb == o_comb)
			return true;

	return false;
}

void MakeCombinations(int n, int r, PermutsVect& output)
{
	output.clear();
	if (n == 0 || r == 0) {
		output.push_back(Permuts());
		return;
	}

	Permuts perm(n);
	for (int i = 0; i < n; ++i)
		perm[i] = i;

	if (n == r) {
		output.push_back(perm);
		return;
	}

	do {
		if (!IsCombDuplicated(r, perm, output)) {
			Permuts comb = perm;
			comb.resize(r);
			std::sort(comb.begin(), comb.end());
			output.push_back(comb);
		}
	} while (std::next_permutation(perm.begin(), perm.end()));
}

void FastMakeCombinations(int n, int r, PermutsVect& output_combs)
{
	output_combs.clear();
	if (n == 0 || r == 0) {
		output_combs.push_back(Permuts());
		return;
	}

	std::vector<int> comb(r);
	for (int i = 0; i < r; ++i)
		comb[i] = i;
	output_combs.push_back(comb);
	if (n == r)
		return;

	while (true) {
		for (int i = comb[r - 1] + 1; i < n; ++i) {
			++comb[r - 1];
			output_combs.push_back(comb);
		}

		int idx = r - 2;
		for (; idx >= 0; --idx) {
			const int MAX_COUNTER = n - r + idx;
			if (comb[idx] < MAX_COUNTER)
				break;
		}

		if (idx < 0)
			break;

		++comb[idx];
		++idx;
		for (; idx < r; ++idx)
			comb[idx] = comb[idx - 1] + 1;
		output_combs.push_back(comb);
	}
}

bool DumpPermutations(const char* file_name, int n, const PermutsVect& perms)
{
	std::ofstream o(file_name);
	if (!o.is_open())
		return false;

	o << n << std::endl << std::endl;
	o << perms.size() << std::endl;
	for (const auto& p : perms) {
		for (int i : p)
			o << i << ' ';
		o << std::endl;
	}
	return true;
}

bool DumpCombinations(const char* file_name, int n, int r, const PermutsVect& combs)
{
	std::ofstream o(file_name);
	if (!o.is_open())
		return false;

	o << n << ' ' << r << std::endl << std::endl;
	o << combs.size() << std::endl;
	for (const auto& c : combs) {
		for (int i : c)
			o << i + 1 << ' ';
		o << std::endl;
	}
	return true;
}

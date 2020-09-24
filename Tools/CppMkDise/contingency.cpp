#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>

#include "asa147.hpp"
#include "libs.h"

//#define READ_K_1
//#define READ_K_2



constexpr int NUM_ALLELES = 3;



typedef std::vector<unsigned char> Chrom;
typedef std::vector<Chrom> Indivs;

typedef std::vector<unsigned char> Labels;



static bool ReadGenos(const char* genos_path, Indivs& indivs)
{
	std::ifstream infile(genos_path);
	if (!infile.is_open()) {
		std::cout << "Could not open genotype `" << genos_path << "' file!" << std::endl;
		return false;
	}

	std::string line;
	int num_indivs, num_loci, num_clusters;
	infile >> line >> num_indivs >> line >> num_loci >> line >> num_clusters;

#if defined FULL_OUTPUT
	std::cout << "NumIndiv:" << num_indivs << std::endl;
	std::cout << "NumLoci :" << num_loci << std::endl;
	std::cout << "NumClus :" << num_clusters << std::endl;
#endif

#if defined READ_K_1 || defined READ_K_2
	indivs.resize(num_indivs);
#else
	indivs.resize(num_indivs * num_clusters);
#endif

	int i = 0;
#if defined READ_K_2
	int num_first = 0;
#endif
	while (std::getline(infile, line)) {
		if (!line.size())
			continue;

#if defined READ_K_2
		if (num_first < num_indivs) {
			++num_first;
			continue;
		}
#endif

		indivs[i].resize(num_loci);
		for (int l = 0; l < num_loci; ++l)
			indivs[i][l] = line[l] - 0x30;

		++i;
#if defined READ_K_1
		if (i >= num_indivs)
			break;
#endif

	}
	return true;
}

static bool ReadLabels(const char* label_path, Labels& lbls, int& num_cases, int& num_controls)
{
	std::ifstream infile(label_path);
	if (!infile.is_open()) {
		std::cout << "Could not open label `" << label_path << "' file!" << std::endl;
		return false;
	}

	std::string line;
	num_cases = num_controls = 0;
	int i = 0;
	while (std::getline(infile, line)) {
#if defined READ_K_2
		if (i < 500) {
			++i;
			continue;
		}
#endif

		const unsigned char LBL = line == "YES" ? 1 : 0;
		if (LBL)
			++num_cases;
		else
			++num_controls;
		lbls.push_back(LBL);

#if defined READ_K_1
		++i;
		if (i >= 500)
			break;
#endif
	}
	return true;
}



double gammainc(double a, double x)
{
	int e;
	const double RES = gammds(x, a, &e);
	if (e) {
#if defined FULL_OUTPUT
		std::cout << "    Error(" << e << ":" <<
			(e == 1 ? "X <= 0 or P <= 0" : "underflow during the computation") <<
			")!    gammds(a:" << a << ", x:" << x << ')' << std::endl;
#endif
		return 1e-10;
	}
	return RES;
}

static double Chi2CDF(double k, double sm)
{
	const double K_2 = k / 2.0;
	const double RES = 1.0 / tgamma(K_2) * gammainc(K_2, sm / 2.0);
	return RES;
}

static double GetPVal(int* cases, int* contr, int num_indivs, int num_cases, int num_controls)
{
	int k = 0;
	double sm = 0.0;
	for (int i = 0; i < NUM_ALLELES; ++i) {
		if (cases[i] < 5 || contr[i] < 5)
			continue;

		++k;
		const double COL_TOT = cases[i] + contr[i];

		const double O_c = cases[i];
		const double E_c = (COL_TOT * num_cases) / num_indivs;
		sm += pow(O_c - E_c, 2) / E_c;

		const double O_u = contr[i];
		const double E_u = (COL_TOT * num_controls) / num_indivs;
		sm += pow(O_u - E_u, 2) / E_u;
	}

	const int DF = k - 1 > 0 ? k - 1 : 1;
	const double RES = 1.0 - Chi2CDF(DF, sm);
	return -log(RES);
}

int main(int argc, char** argv)
{
#if defined FULL_OUTPUT
	std::cout << GetCTime() << "Start contingency . . ." << std::endl;
#endif

	Indivs indivs;
	if (!ReadGenos("cpp_dise_genos_K2_L10000_D5_N500.str-faststr.txt", indivs))
		return 1;

	const int NUM_INDIVS = static_cast<int>(indivs.size());
	const int NUM_LOCI = static_cast<int>(indivs[0].size());

	Labels labels;
	int num_cases, num_controls;
	if (!ReadLabels("cpp_dise_label_K2_L10000_D5_N500.txt", labels, num_cases, num_controls))
		return 2;


	int cases[NUM_ALLELES];
	int contr[NUM_ALLELES];
	for (int l = 0; l < NUM_LOCI; ++l) {
		memset(&cases, 0, sizeof(cases[0]) * NUM_ALLELES);
		memset(&contr, 0, sizeof(cases[0]) * NUM_ALLELES);
		for (int i = 0; i < NUM_INDIVS; ++i) {
			const int IDX = static_cast<int>(indivs[i][l]);
			if (labels[i])
				++cases[IDX];
			else
				++contr[IDX];
		}
		double p_val = GetPVal(cases, contr, NUM_INDIVS, num_cases, num_controls);
#if defined FULL_OUTPUT
		std::cout << l <<
			"     C ->   " << cases[0] << ' ' << cases[1] << ' ' << cases[2] <<
			"     U ->   " << contr[0] << ' ' << contr[1] << ' ' << contr[2] <<
			"     p-val " << p_val << std::endl;
#else
		std::cout << ((l + 1) * 0.1) << '\t' << p_val << std::endl;
#endif
	}
	return 0;
}


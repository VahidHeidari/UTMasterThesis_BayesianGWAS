#include <sstream>

#include "bit_genos_matrix.h"
#include "dists.h"
#include "logger.h"
#include "params.h"



static void TestFastMakeCombinations()
{
	const int N = 10;
	for (int R = 0; R <= N; ++R) {
		PermutsVect combs1;
		std::stringstream ss;
		ss << BASE_PATH << "C-C-" << N << '-' << R << ".txt";
		FastMakeCombinations(N, R, combs1);
		DumpCombinations(ss.str().c_str(), N, R, combs1);
		logger << "N:" << N << "   R:" << R << std::endl;
		logger << "End COUNTER    : " << Time << std::endl;

		PermutsVect combs2;
		ss.str("");
		ss << BASE_PATH << "C-E-" << N << '-' << R << ".txt";
		MakeCombinations(N, R, combs2); DumpCombinations(ss.str().c_str(), N, R, combs2);
		logger << "End EXHUASTIVE : " << Time << std::endl;

		bool is_ok = true;
		if (combs1.size() != combs2.size()) {
			logger << warning << "Size is not equal! counter(" << combs1.size() << ")    exhuastive(" << combs2.size() << ")" << std::endl;
			is_ok = false;
		}
		for (unsigned i = 0; i < combs1.size(); ++i)
			for (unsigned j = 0; j < combs1[i].size(); ++j)
				if (combs1[i][j] != combs2[i][j]) {
					logger << warning << "combs1[" << i << "][" << j << "](" << combs1[i][j]
							<< ") != combs2[" << i << "][" << j << "](" << combs2[i][j] << ")" << std::endl;
					is_ok = false;
					goto log_lbl;
				}
log_lbl:
		if (is_ok)
			logger << "every thing is OK" << std::endl;
		else
			return;

		logger << std::endl;
	}
}



static void TestGenosMatrix()
{
	// First test
	static const unsigned int TEST_GENOS_1[5 * 5] =
	{
		2, 1, 0, 2, 2,
		2, 2, 0, 0, 1,
		1, 1, 1, 2, 1,
		0, 2, 2, 2, 0,
		1, 0, 0, 1, 2,
	};

	bool is_ok = true;
	BitGenosMatrix genos(5, 5, 1);
#define SET_GENO(I, L, G)	genos.SetGeno(I, L, G)
	SET_GENO(0, 0, 2); SET_GENO(0, 1, 1); SET_GENO(0, 2, 0); SET_GENO(0, 3, 2); SET_GENO(0, 4, 2);
	SET_GENO(1, 0, 2); SET_GENO(1, 1, 2); SET_GENO(1, 2, 0); SET_GENO(1, 3, 0); SET_GENO(1, 4, 1);
	SET_GENO(2, 0, 1); SET_GENO(2, 1, 1); SET_GENO(2, 2, 1); SET_GENO(2, 3, 2); SET_GENO(2, 4, 1);
	SET_GENO(3, 0, 0); SET_GENO(3, 1, 2); SET_GENO(3, 2, 2); SET_GENO(3, 3, 2); SET_GENO(3, 4, 0);
	SET_GENO(4, 0, 1); SET_GENO(4, 1, 0); SET_GENO(4, 2, 0); SET_GENO(4, 3, 1); SET_GENO(4, 4, 2);
#undef SET_GENO
	for (int i = 0; i < genos.GetNumIndivs(); ++i) {
		for (int l = 0; l < genos.GetNumLoci(); ++l) {
			const int g = genos.GetGeno(i, l);
			if (g != TEST_GENOS_1[i * genos.GetNumLoci() + l])
				is_ok = false;
			logger << g << ' ';
		}
		logger << std::endl;
	}

	if (is_ok)
		logger << "Every thing is OK!" << std::endl;
	else
		logger << "Some thing is going worng!" << std::endl;


	// Second test
	static const unsigned int TEST_GENOS_2[17 * 3] =
	{
		2, 1, 0,
		2, 2, 0,
		1, 1, 1,
		2, 1, 0,
		0, 0, 0,
		0, 1, 2,
		2, 1, 0,
		1, 2, 1,
		2, 1, 2,
		0, 1, 0,
		1, 0, 1,
		1, 2, 1,
		2, 0, 2,
		0, 2, 1,
		2, 2, 2,
		0, 2, 1,
		1, 1, 0,
	};
	BitGenosMatrix genos2(17, 3, 1);
#define SET_GENO(I, L, G)	genos2.SetGeno(I, L, G)
	SET_GENO( 0, 0, 2); SET_GENO( 0, 1, 1); SET_GENO( 0, 2, 0);
	SET_GENO( 1, 0, 2); SET_GENO( 1, 1, 2); SET_GENO( 1, 2, 0);
	SET_GENO( 2, 0, 1); SET_GENO( 2, 1, 1); SET_GENO( 2, 2, 1);
	SET_GENO( 3, 0, 2); SET_GENO( 3, 1, 1); SET_GENO( 3, 2, 0);
	SET_GENO( 4, 0, 0); SET_GENO( 4, 1, 0); SET_GENO( 4, 2, 0);
	SET_GENO( 5, 0, 0); SET_GENO( 5, 1, 1); SET_GENO( 5, 2, 2);
	SET_GENO( 6, 0, 2); SET_GENO( 6, 1, 1); SET_GENO( 6, 2, 0);
	SET_GENO( 7, 0, 1); SET_GENO( 7, 1, 2); SET_GENO( 7, 2, 1);
	SET_GENO( 8, 0, 2); SET_GENO( 8, 1, 1); SET_GENO( 8, 2, 2);
	SET_GENO( 9, 0, 0); SET_GENO( 9, 1, 1); SET_GENO( 9, 2, 0);
	SET_GENO(10, 0, 1); SET_GENO(10, 1, 0); SET_GENO(10, 2, 1);
	SET_GENO(11, 0, 1); SET_GENO(11, 1, 2); SET_GENO(11, 2, 1);
	SET_GENO(12, 0, 2); SET_GENO(12, 1, 0); SET_GENO(12, 2, 2);
	SET_GENO(13, 0, 0); SET_GENO(13, 1, 2); SET_GENO(13, 2, 1);
	SET_GENO(14, 0, 2); SET_GENO(14, 1, 2); SET_GENO(14, 2, 2);
	SET_GENO(15, 0, 0); SET_GENO(15, 1, 2); SET_GENO(15, 2, 1);
	SET_GENO(16, 0, 1); SET_GENO(16, 1, 1); SET_GENO(16, 2, 0);
#undef SET_GENO
	is_ok = true;
	for (int i = 0; i < genos2.GetNumIndivs(); ++i) {
		for (int l = 0; l < genos2.GetNumLoci(); ++l) {
			const int G = genos2.GetGeno(i, l);
			const int EXP_G = TEST_GENOS_2[i * genos2.GetNumLoci() + l];
			if (G != EXP_G)
				is_ok = false;
			logger << G << ' ';
		}
		logger << std::endl;
	}
	if (is_ok)
		logger << "Every thing is OK!" << std::endl;
	else
		logger << "Some thing is going worng!" << std::endl;
}



int main()
{
	logger << std::endl << std::endl;
	logger << "----- UNIT TEST -----" << std::endl;
	logger << "Start : " << Time << std::endl;

	//TestFastMakeCombinations();
	TestGenosMatrix();

	logger << "End : " << Time << std::endl << std::endl;
	logger << " [OK] TEST PASSED!" << std::endl << std::endl;
	return 0;
}

#include "libs.h"

//#define ALL_MODELS 1

#if defined ALL_MODELS
static std::vector<DisModel> models = {
	{
		{{220,  0}}, {{2392, 1}}, {{5505, 2}},		// Single locus
		{{6183, 0}, {8139, 2}},						// Two way interaction
		{{539,  2}, {3242, 0}},
		{{1724, 0}, {2620, 0}, {6005, 1}}			// Three way interaction
	}, {
		{{78,   1}}, {{7389, 1}}, {{4216, 1}},		// Single locus
		{{5435, 0}, {9491, 1}},						// Two way interaction
		{{2969, 2}, {8679, 1}},
		{{578,  1}, {3786, 2}, {9970, 1}}			// Three way interaction
	},
};

static Params params = {
	{ -0.1,    0.5, 0.5, 0.5,    1.0, 1.0,    1.5 },
	{ -0.1,    0.5, 0.5, 0.5,    1.0, 1.0,    1.5 },
};
#else
static std::vector<DisModel> models = {
	{
		{{1724, 0}, {2620, 0}, {6005, 1}}			// Three way interaction
	}, {
		{{578,  1}, {3786, 2}, {9970, 1}}			// Three way interaction
	},
};

static Params params = {
	{ -0.1, 1.5 },
	{ -0.1, 1.5 },
};
#endif


int main(int argc, char** argv)
{
	if (argc < 3) {
		std::cout << " Usage: mk_dise NUM_INDIVS DIFF_PERCENT" << std::endl;
		return 1;
	}

	NUM_LOCI = 10000;
	NUM_CLUSTERS = 2;
	NUM_INDIVS = std::stoi(argv[1]);
	const int DIFF_PERCENT = std::stoi(argv[2]);
	const int DIFFS_LOCI = static_cast<int>(static_cast<double>(DIFF_PERCENT) / 100.0 * NUM_LOCI + 0.5);
	NUM_DIFFS = (DIFFS_LOCI <= 0) ? 1 : DIFFS_LOCI;

	FMT =
		"_K" + std::to_string(NUM_CLUSTERS) +
		"_L" + std::to_string(NUM_LOCI) +
		"_D" + std::to_string(DIFF_PERCENT) +
		"_N" + std::to_string(NUM_INDIVS);

	FREQS_PATH = std::string("cpp_dise_freqs") + FMT + ".txt";
	GENOS_PATH = std::string("cpp_dise_genos") + FMT + ".str";
	LABEL_PATH = std::string("cpp_dise_label") + FMT + ".txt";

	std::cout << "  MK DISE" << std::endl;
	std::cout << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	std::cout << "NUM_DIFFS:    " << NUM_DIFFS << std::endl;
	std::cout << "NUM_INDIVS:   " << NUM_INDIVS << std::endl;

	std::cout << GetCTime() << "START" << std::endl;
	FreqsVect freqs;
	MakeDiseaseFreqs(freqs, models);
	MakeDiseGenos(freqs, models, params);
	std::cout << GetCTime() << "END" << std::endl;
	return 0;
}


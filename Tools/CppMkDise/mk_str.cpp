#include "libs.h"

int main(int argc, char** argv)
{
	if (argc < 5) {
		std::cout << " Usage: mk_str NUM_CLUSTERS NUM_LOCI NUM_INDIVS DIFF_PERCENT" << std::endl;
		return 1;
	}

	NUM_CLUSTERS = std::stoi(argv[1]);
	NUM_LOCI = std::stoi(argv[2]);
	NUM_INDIVS = std::stoi(argv[3]);
	const int DIFF_PERCENT = std::stoi(argv[4]);
	const int DIFFS_LOCI = static_cast<int>(static_cast<double>(DIFF_PERCENT) / 100.0 * NUM_LOCI + 0.5);
	NUM_DIFFS = (DIFFS_LOCI <= 0) ? 1 : DIFFS_LOCI;

	FMT =
		"_K" + std::to_string(NUM_CLUSTERS) +
		"_L" + std::to_string(NUM_LOCI) +
		"_D" + std::to_string(DIFF_PERCENT) +
		"_N" + std::to_string(NUM_INDIVS);

	FREQS_PATH = std::string("cpp_freqs") + FMT + ".txt";
	GENOS_PATH = std::string("cpp_genos") + FMT + ".str";

	std::cout << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	std::cout << "NUM_LOCI:     " << NUM_LOCI << std::endl;
	std::cout << "NUM_DIFFS:    " << NUM_DIFFS << std::endl;
	std::cout << "NUM_INDIVS:   " << NUM_INDIVS << std::endl;

	std::cout << GetCTime() << "START" << std::endl;
	FreqsVect freqs;
	MakeFreqs(freqs);
	MakeGenos(freqs);
	std::cout << GetCTime() << "END" << std::endl;
	return 0;
}


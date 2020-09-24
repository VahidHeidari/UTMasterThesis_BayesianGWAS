#include <string>

#include "logger.h"
#include "SAM-to-genotypes.h"

static std::string GetFileName(const std::string& file_path)
{
	auto last_slash_pos = file_path.rfind('\\');
	return last_slash_pos == std::string::npos ? file_path : file_path.substr(last_slash_pos + 1);
}

int main(int argc, char** argv)
{
	if (argc < 3) {
		logger << Time << "Fasta file and SAM file name are required!" << std::endl;
		return 1;
	}

	// Process import command line options.
	std::string fasta_path = { argv[1] };
	std::string sam_path = { argv[2] };
	std::string out_path = GetFileName(sam_path) + ".gtp";

	// Process optional command line options.
	if (argc > 3)
		out_path = std::string(argv[3]) + '\\' + out_path;

#ifndef USE_LAZY_GET_NUCLEOTIDE
	FastaFile fasta;
	if (!fasta.ReadFile(fasta_path)) {
		logger << Time << "Could not read Fasta file!" << std::endl;
		return 2;
	}
	logger << Time << "Fasta size : " << fasta.GetBytes() << " bytes (" << GetHumanSize(fasta.GetBytes()) << ")" << std::endl;
#else
	FastaFile fasta(fasta_path);
#endif

	WriteToGenotype(sam_path, out_path, fasta);
	return 0;
}

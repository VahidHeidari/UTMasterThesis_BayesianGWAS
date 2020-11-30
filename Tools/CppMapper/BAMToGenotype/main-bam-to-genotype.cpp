#include <string>

#include "BAM-file.h"
#include "Fasta-file.h"
#include "logger.h"
#include "SAM-to-genotypes.h"

static std::string GetFileName(const std::string& file_path)
{
	auto last_slash_pos = file_path.rfind('\\');
	return last_slash_pos == std::string::npos ? file_path : file_path.substr(last_slash_pos + 1);
}

int main(int argc, char** argv)
{
	logger << std::endl << std::endl << Time << "Start of BAM2Genotype" << std::endl;
	if (argc < 3) {
		logger << Time << "Fasta file and BAM file name are required!" << std::endl;
		return 1;
	}

	// Process import command line options.
	std::string fasta_path = { argv[1] };
	std::string bam_path = { argv[2] };
	std::string out_path = GetFileName(bam_path) + ".gtp";

	// Process optional command line options.
	if (argc > 3)
		out_path = std::string(argv[3]) + '\\' + out_path;

#ifndef USE_LAZY_GET_NUCLEOTIDE
	logger << Time << "Loading full Fasta file . . ." << std::endl;
	FastaFile fasta;
	if (!fasta.ReadFile(fasta_path)) {
		logger << Time << "Could not read Fasta file!" << std::endl;
		return 2;
	}
	logger << Time << "Fasta size : " << fasta.GetBytes() <<
		" bytes (" << GetHumanSize(fasta.GetBytes()) << ")" << std::endl;
#else
	logger << Time << "Loading lazy Fasta file . . ." << std::endl;
	FastaFile fasta(fasta_path);
#endif
	logger << Time << "FASTA input path  -> " << fasta_path << std::endl;

	WriteBAMToGenotype(bam_path, out_path, fasta);
	logger << Time << "END of execution!" << std::endl << std::endl;
	return 0;
}


#include "SAM-to-genotypes.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "logger.h"
#include "SAM-file.h"



static long long Pow(int a, int b)
{
	long long res = 1;
	for (int i = 0; i < b; ++i)
		res *= a;
	return res;
}

std::string GetHumanSize(long long sz)
{
	static const char* SZ_STR[] = { "Bytes", "KB", "MB", "GB", "TB" };
	int i = 0;
	while (sz / Pow(1000, i + 1))
		++i;
	double szf = static_cast<double>(sz) / static_cast<double>(Pow(1024, i));
	std::ostringstream ss;
	ss << std::fixed << std::setprecision(1) << szf << ' ' << SZ_STR[i];
	return ss.str();
}

static void DumpGenotypeList(GenotypeList& genotype,
	std::ostream& genotype_file, const std::string& cur_rname)
{
	if (cur_rname == "MT")
		return;

	logger << Time << "Dump reference `" << cur_rname << "' to output" << std::endl;
	if (IS_BINARY_OUTPUT)
		genotype_file.write(cur_rname.c_str(), cur_rname.size() + 1);
	else
		genotype_file << ">>> " << cur_rname << std::endl;
	genotype.Dump(genotype_file);
	logger << Time << "Dump completed!" << std::endl;
}

void WriteToGenotype(const std::string& sam_path, const std::string& out_path,
#ifndef USE_LAZY_GET_NUCLEOTIDE
	const FastaFile& fasta
#else
	FastaFile& fasta
#endif
)
{
	logger << std::endl << std::endl << std::endl << Time << std::endl;
	logger << "---------- START ----------" << std::endl;
	logger << Time << "SAM input path    -> " << sam_path << std::endl;
	logger << Time << "Genotype out path -> " << out_path << std::endl;
	if (IS_BINARY_OUTPUT || IS_COMPRESSED_OUTPUT)
		logger << Time << "Output is binary and " << (IS_COMPRESSED_OUTPUT ? "" : "not ") << "compressed" << std::endl;
	else
		logger << Time << "Output is text" << std::endl;

	GenotypeList genotype;
	std::ifstream sam_file(sam_path);
	if (!sam_file.is_open()) {
		logger << Time << "Could not open input file!" << std::endl;
		return;
	}

	std::ofstream genotype_file(out_path, OUT_FILE_FLAGS);
	if (!genotype_file.is_open()) {
		logger << Time << "Could not open output file!" << std::endl;
		return;
	}

	long long num_sam_recs = 0;
	long long tot_loci = 0;
	CIGAROps ops;
	std::string line;
	std::string cur_rname = "";
	while (std::getline(sam_file, line)) {
		if (line[0] == '@')
			continue;

		++num_sam_recs;
		if (num_sam_recs % NUM_SAM_REC_PRINT == 0)
			logger << Time << "Number of SAM records: " << num_sam_recs
			<< "    mem size: " << GetHumanSize(genotype.GetMemorySize())
			<< "    num loci: " << genotype.GetNumLoci() << std::endl;

		// Split SAM record.
		const auto rec = Split(line);
		if (HasFlag(rec, FLAGS_UNMAPPED))
			continue;

		const std::string& rname = GetRName(rec);
		if (rname != cur_rname) {
			if (cur_rname.size()) {
				logger << Time << "NEW CHROMOSOME `" << rname << "'!" << std::endl;
				if (rname == "MT") {
					logger << Time << "Skip MT chromosome!" << std::endl;
					cur_rname = rname;
					break;
				}

				logger << Time << "Number of SAM records: " << num_sam_recs
					<< "    mem size: " << GetHumanSize(genotype.GetMemorySize())
					<< "    num loci: " << genotype.GetNumLoci() << std::endl;

				// Dump to file and reset for the next chromosome.
				DumpGenotypeList(genotype, genotype_file, cur_rname);
				tot_loci += genotype.GetNumLoci();
				genotype.Reset();
			}
			cur_rname = rname;
		}

		if (GetCIGAR(rec) == "*")
			continue;

		// Find genotype list position.
		const std::string& seq = GetSeq(rec);
		int seq_idx = 0, start_sam_pos, sam_pos;
		start_sam_pos = sam_pos = GetPos(rec);
		SplitCIGAR(rec, ops);
		for (size_t i = 0; i < ops.ops.size(); ++i) {
			if (ops.ops[i] == 'H' || ops.ops[i] == 'P')
				continue;

			if (ops.ops[i] == 'S') {
				seq_idx += ops.ops_len[i];
				continue;
			}

			if (ops.ops[i] == 'D' || ops.ops[i] == 'N') {
				sam_pos += ops.ops_len[i];
				continue;
			}

			// Insertion sequence.
			const int SUB_SEQ_LEN = ops.ops_len[i];
			if (ops.ops[i] == 'I') {
				const int j_start = (i - 1 >= 0 && ops.ops[i - 1] == 'P') ? ops.ops_len[i - 1] : 0;
				for (int j = j_start; j < SUB_SEQ_LEN; ++j)
					genotype.AddInsert(sam_pos, j, seq[seq_idx + j]);
				seq_idx += ops.ops_len[i];
				continue;
			}

			// Matched sequences.
			for (int j = 0; j < SUB_SEQ_LEN; ++j)
				genotype.Add(sam_pos++, seq[seq_idx + j], ops.ops[i], cur_rname, fasta);
			seq_idx += ops.ops_len[i];
		}
	}

	// Dump last chromosome.
	DumpGenotypeList(genotype, genotype_file, cur_rname);
	tot_loci += genotype.GetNumLoci();

	// Close files.
	sam_file.close();
	genotype_file.close();

	logger << Time << "Number of SAM records: " << num_sam_recs << std::endl;
	logger << Time << "Number of loci       : " << tot_loci << std::endl;
	logger << Time << "Done!" << std::endl;
}

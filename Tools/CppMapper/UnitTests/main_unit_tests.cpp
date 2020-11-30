#include <fstream>
#include <iostream>
#include <vector>

#include "SAM-to-genotypes.h"


constexpr int GENOTYPE_MASK = 0x03;



#define CHECK_ITEMS_VECT(KEY, IDX, ...)	do {	\
	std::vector<int> items_vect { __VA_ARGS__ };		\
	const int idx = BinSearch(KEY, items_vect);			\
	if (idx == IDX)										\
		std::cout << "OK!" << std::endl;				\
	else												\
		std::cout << "FAILED! EXP:" << IDX				\
			<< " idx:" << idx << std::endl				\
			<< "LINE:" << __LINE__ << ' '				\
			<< "FILE:" << __FILE__ << std::endl;		\
} while (false)



static int BinSearch(int key, const std::vector<int>& items_vect)
{
	if (!items_vect.size())
		return 0;

	if (items_vect.size() == 1)
		return (items_vect[0] >= key) ? 0 : 1;

	int left = 0, right = static_cast<int>(items_vect.size()) - 1;
	while (left < right) {
		if (items_vect[right] == key)
			return right;

		if (items_vect[left] == key)
			return left;

		const int mid = left + (right - left) / 2;
		if (items_vect[mid] == key)
			return mid;

		if (items_vect[mid] < key)
			left = mid + 1;
		else
			right = mid - 1;
	}

	return right >= 0 && items_vect[right] > key ? right : right + 1;
}

static void PrintItemsVect(const std::vector<int>& items_vect)
{
	for (auto& r : items_vect)
		std::cout << r << ' ';
	std::cout << std::endl;
}

static void TestBinarySearch()
{
	// Single item test
	CHECK_ITEMS_VECT(1, 0, 1);
	CHECK_ITEMS_VECT(0, 0, 1);
	CHECK_ITEMS_VECT(2, 1, 1);

	// Two items test
	CHECK_ITEMS_VECT(0, 0, 1, 2);
	CHECK_ITEMS_VECT(1, 0, 1, 2);
	CHECK_ITEMS_VECT(2, 1, 1, 2);
	CHECK_ITEMS_VECT(3, 2, 1, 2);

	// Three items test
	CHECK_ITEMS_VECT(0, 0, 1, 2, 3);
	CHECK_ITEMS_VECT(1, 0, 1, 2, 3);
	CHECK_ITEMS_VECT(2, 1, 1, 2, 3);
	CHECK_ITEMS_VECT(3, 2, 1, 2, 3);
	CHECK_ITEMS_VECT(4, 3, 1, 2, 3);
	CHECK_ITEMS_VECT(9, 3, 1, 2, 3);

	// Random test
	CHECK_ITEMS_VECT(2, 0, 3);
	CHECK_ITEMS_VECT(2, 1, 1, 3);
	CHECK_ITEMS_VECT(2, 1, 1);
	CHECK_ITEMS_VECT(4, 2, 1, 3);

	// One hole in middel test
	CHECK_ITEMS_VECT(0, 0, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(1, 0, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(2, 1, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(3, 1, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(4, 2, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(5, 2, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(6, 3, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(7, 3, 1, 3, 5, 7);
	CHECK_ITEMS_VECT(8, 4, 1, 3, 5, 7);
}



static void ReadCompressedGenotype(std::string gtp_file_path)
{
	std::ifstream gtp_file(gtp_file_path);
	if (!gtp_file.is_open()) {
		logger << Time << "Could not open gtp file!" << std::endl;
		return;
	}

	// Read reference name.
	std::string rname;
	while (true) {
		char ch;
		gtp_file >> ch;
		if (ch == '\0')
			break;

		rname += ch;
	}
	std::cout << "UnitTests ReadCompressedGenotype() -> ref name is: " << rname << std::endl;

	uint32_t num_genotypes, num_inserts;
	gtp_file.read(reinterpret_cast<char*>(&num_genotypes), sizeof(num_genotypes));
	gtp_file.read(reinterpret_cast<char*>(&num_inserts), sizeof(num_inserts));
	std::cout << "num geno: " << num_genotypes << std::endl;
	std::cout << "num insr: " << num_inserts << std::endl;

	for (unsigned i = 0; i < num_genotypes; ++i) {
		uint32_t pos, len;
		gtp_file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
		gtp_file.read(reinterpret_cast<char*>(&len), sizeof(len));
		std::cout << "pos:" << pos << "     len:" << len << std::endl;
		int bit_num = 0;
		uint32_t buff;
		gtp_file.read(reinterpret_cast<char*>(&buff), sizeof(buff));
		for (unsigned j = 0; j < len; ++j) {
			if (bit_num == BUFFER_BITS) {
				bit_num = 0;
				gtp_file.read(reinterpret_cast<char*>(&buff), sizeof(buff));
			}
			std::cout << ((buff >> bit_num) & GENOTYPE_MASK) << ' ';
			bit_num += BITS_PER_GENOTYPE;
		}
		std::cout << std::endl;
		i += len;
	}
	std::cout << std::endl;
}

static void TestSAMToGenotypeCompressed()
{
#ifdef USE_LAZY_GET_NUCLEOTIDE
	FastaFile fasta("D:\\Datasets\\Programs\\Mapper\\SAMSpecSamples\\SAM_BAM_spec.fasta.bin");
#else
	FastaFile fasta;
	if (!fasta.ReadFile("D:\\Datasets\\Programs\\Mapper\\SAMSpecSamples\\SAM_BAM_spec.fasta.bin")) {
		std::cout << "Could not read test Fasta file!" << std::endl;
		return;
	}
#endif
	const char REF[] = "AGCATGTTAGATAAGATAGCTGTGCTAGTAGGCAGTCAGCGCCAT";
	bool is_ok = true;
	for (int i = 0; i < 45; ++i) {
		const char ch = fasta.GetNucleotide("ref", i + 1);
		std::cout << ch;
		if (ch != REF[i])
			is_ok = false;
	}
	std::cout << std::endl;
	std::cout << "Reference is " << (is_ok ? "OK!" : "Not OK!") << std::endl;

	WriteToGenotype("D:\\Datasets\\Programs\\Mapper\\SAMSpecSamples\\SAM_BAM_spec.sam",
		"D:\\Datasets\\Programs\\CppMapper\\SAM_BAM_spec.gtp", fasta);

	ReadCompressedGenotype("D:\\Datasets\\Programs\\CppMapper\\SAM_BAM_spec.gtp");
}

static void DumpFastaFile()
{
	//const std::string input_path = "D:\\Datasets\\Programs\\Mapper\\SAMSpecSamples\\SAM_BAM_spec.fasta";
	const std::string input_path = "D:\\Datasets\\Homo_sapiens.GRCh38.dna.primary_assembly.fa";
	const std::string output_path = input_path + ".bin";

	FastaFile fasta;
	if (!fasta.ReadFile(input_path)) {
		std::cout << "Could not open fasta file for dumping!" << std::endl;
		return;
	}

	if (!fasta.DumpBinary(output_path))
		std::cout << "Could not dump fasta file!" << std::endl;
}



inline static std::string ReadRefName(std::ifstream& f)
{
	std::string rname;
	char ch;
	while (true) {
		f.read(&ch, 1);
		if (ch == '\0')
			break;

		rname += ch;
	}
	return rname;
}

inline static std::ifstream::pos_type GetFileSize(const std::string& file_path)
{
	std::ifstream f(file_path, std::ios_base::in | std::ios_base::out | std::ios_base::binary);
	if (!f.is_open())
		return std::fstream::pos_type(-1);

	f.seekg(0, std::ios_base::end);
	return f.tellg();
}

static void ConvertCompressedGTPToTextGTP(const std::string& comp_gtp_path,
		const std::string& out_path)
{
	std::ifstream gtp_file(comp_gtp_path, std::ios_base::in | std::ios_base::binary);
	if (!gtp_file.is_open()) {
		logger << Time << "Could not open gtp file!" << std::endl;
		return;
	}

	std::ofstream gtp_out_file(out_path);
	if (!gtp_out_file.is_open()) {
		logger << Time << "Could not open gtp output file!" << std::endl;
		return;
	}

	std::ifstream::pos_type gtp_file_size = GetFileSize(comp_gtp_path);
	std::string rname;
	while (gtp_file.tellg() < gtp_file_size) {
		rname = ReadRefName(gtp_file);
		logger << Time << "UnitTests ReadCompressedGenotype() -> ref name is: " << rname
			<< "    " << gtp_file.tellg() << std::endl;
		gtp_out_file << ">>> " << rname << std::endl;

		uint32_t num_genotypes, num_inserts;
		gtp_file.read(reinterpret_cast<char*>(&num_genotypes), sizeof(num_genotypes));
		gtp_file.read(reinterpret_cast<char*>(&num_inserts), sizeof(num_inserts));
		logger << Time  << "num geno: " << num_genotypes << std::endl;
		logger << Time  << "num insr: " << num_inserts << std::endl;

		for (unsigned i = 0; i < num_genotypes; ++i) {
			uint32_t pos, len;
			gtp_file.read(reinterpret_cast<char*>(&pos), sizeof(pos));
			gtp_file.read(reinterpret_cast<char*>(&len), sizeof(len));
			//logger << Time  << "pos:" << pos << "     len:" << len << std::endl;
			int bit_num = 0;
			uint32_t buff;
			gtp_file.read(reinterpret_cast<char*>(&buff), sizeof(buff));
			for (unsigned j = 0; j < len; ++j) {
				if (bit_num == BUFFER_BITS) {
					bit_num = 0;
					gtp_file.read(reinterpret_cast<char*>(&buff), sizeof(buff));
				}
				//logger << ' ' << pos + j << ' ' << "ACTG"[((buff >> bit_num) & GENOTYPE_MASK)] << std::endl;
				gtp_out_file <<  ' ' << pos + j << ' ' << "ACTG"[((buff >> bit_num) & GENOTYPE_MASK)] << std::endl;
				bit_num += BITS_PER_GENOTYPE;
			}
			//logger << std::endl;
			i += len;
		}
		//logger << std::endl;
		logger << Time << "Try break" << std::endl;
		break;
	}
}



int main()
{
	//TestBinarySearch();
	//TestSAMToGenotypeCompressed();
	//DumpFastaFile();
	ConvertCompressedGTPToTextGTP(
			"D:\\Datasets\\GSE68086\\CRC-sorted\\SRR1982618-sorted.bam.gtp",
			"D:\\Datasets\\GSE68086\\CRC-sorted\\SRR1982618-sorted.bam.gtp.txt");
	return 0;
}

#ifndef FASTA_FILE_H_
#define FASTA_FILE_H_

#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "logger.h"

#define USE_LAZY_GET_NUCLEOTIDE		1

/// Forward decleartion
std::string GetHumanSize(long long sz);

class FastaFile
{
	typedef std::vector<unsigned char> Chromosome;
	typedef std::map<std::string, Chromosome*> ChromeMap;

	enum { BINARY_FILE_SIGNATURE = 0xFFEECCAA };

public:
	FastaFile() {}
	FastaFile(std::string fasta_path) : fasta_path(fasta_path) {}
	~FastaFile()
	{
		for (auto& itr : chromosomes)
			delete itr.second;
	}

	bool ReadFile(std::string fasta_path)
	{
		if (IsBinaryFile(fasta_path))
			return ReadBinaryFile(fasta_path);

		return ReadTextFile(fasta_path);
	}

#ifndef USE_LAZY_GET_NUCLEOTIDE
	char GetNucleotide(const std::string& ref_name, int pos) const
	{
		ChromeMap::const_iterator& itr = chromosomes.find(ref_name);
		if (itr == chromosomes.end()) {
			logger << Time << "WARNING: Could not find `" << ref_name << "' chromosome!" << std::endl;
			return '\0';
		}

		const auto& vect = *(itr->second);
		const char ref = GetRefChar(vect, ref_name, pos);
		return ref;
	}
#else
	char GetNucleotide(const std::string& ref_name, int pos)
	{
		if (cur_name == ref_name)
			return GetRefChar(lazy_chromosome, ref_name, pos);

		logger << Time << "FastaFile: Load new chromosome `" << ref_name << "'" << std::endl;
		lazy_chromosome.clear();
		std::ifstream fasta_file(fasta_path);

		// Signature
		uint32_t sig;
		fasta_file.read(reinterpret_cast<char*>(&sig), sizeof(sig));
		if (sig != BINARY_FILE_SIGNATURE) {
			logger << Time << "FastaFile: The input file is not binary fasta file!" << std::endl;
			return '\0';
		}

		// Number of chromosomes
		uint32_t num_chromosomes;
		fasta_file.read(reinterpret_cast<char*>(&num_chromosomes), sizeof(num_chromosomes));

		// Finding reference name
		for (unsigned i = 0; i < num_chromosomes; ++i) {
			std::string tmp_ref_name = ReadRefName(fasta_file);
			uint32_t num_bytes;
			fasta_file.read(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
			if (tmp_ref_name != ref_name) {
				fasta_file.seekg(num_bytes, std::ios_base::cur);
				continue;
			}

			lazy_chromosome.resize(num_bytes);
			fasta_file.read(reinterpret_cast<char*>(&lazy_chromosome[0]), num_bytes);
			cur_name = ref_name;
			logger << Time << "FastaFile: Chromosome size: " << num_bytes << " bytes ("
				<< GetHumanSize(num_bytes) << ")" << std::endl;
			return GetRefChar(lazy_chromosome, ref_name, pos);
		}
		logger << Time << "FastaFile: lazy reading could not find chromosome!" << std::endl;
		return '\0';
	}
#endif

	bool DumpBinary(std::string out_path) const
	{
		logger << Time << "FastaFile: Dump content of fasta file to binary file . . ." << std::endl;
		std::ofstream out_file(out_path);
		if (!out_file.is_open()) {
			logger << Time << "FastaFile: Could not open output file! Dumping failed." << std::endl;
			return false;
		}

		// Signature
		const uint32_t SIG = BINARY_FILE_SIGNATURE;
		out_file.write(reinterpret_cast<const char*>(&SIG), sizeof(SIG));

		// Number of chromosomes
		const uint32_t NUM_CHROMOSOMES = static_cast<int>(chromosomes.size());
		out_file.write(reinterpret_cast<const char*>(&NUM_CHROMOSOMES), sizeof(NUM_CHROMOSOMES));

		// Chromosomes data
		for (const auto& itr : chromosomes) {
			const auto& ref_name = itr.first;
			out_file.write(ref_name.c_str(), ref_name.size() + 1);
			const uint32_t num_bytes = static_cast<int>(itr.second->size());
			out_file.write(reinterpret_cast<const char*>(&num_bytes), sizeof(num_bytes));
			const Chromosome& vect = itr.second[0];
			out_file.write(reinterpret_cast<const char*>(&vect[0]), num_bytes);
		}
		logger << Time << "FastaFile: Dump is done!" << std::endl;
		return true;
	}

	int GetBytes() const
	{
		int bytes = 0;
		for (const auto& itr : chromosomes)
			bytes += static_cast<int>(itr.second->size());
		return bytes;
	}

private:
	bool IsBinaryFile(const std::string& fasta_path) const
	{
		std::ifstream test_file(fasta_path);
		if (!test_file.is_open())
			return false;

		uint32_t sig;
		test_file.read(reinterpret_cast<char*>(&sig), sizeof(sig));
		return sig == BINARY_FILE_SIGNATURE;
	}

	std::string GetRefName(const std::string& line) const
	{
		unsigned i = 1;
		while (i < line.size()) {
			if (line[i] == ' ' || line[i] == '\t' || line[i] == '\r' || line[i] == '\n')
				break;

			++i;
		}

		std::string rname = line.substr(1, i - 1);
		return rname;
	}

	bool ReadTextFile(const std::string& fasta_path)
	{
		logger << Time << "FastaFile: Reading input -> " << fasta_path << std::endl;
		std::ifstream fasta_file(fasta_path);
		if (!fasta_file.is_open()) {
			logger << Time << "FastaFile: Could not reading text input file!" << std::endl;
			return false;
		}

		std::string cur_rname;
		Chromosome* refs = new Chromosome;
		std::string line;
		while (std::getline(fasta_file, line)) {
			if (line[0] == '>') {
				std::string rname = GetRefName(line);
				if (rname == "MT") {
					logger << Time << "FastaFile: Skip `MT' chromsomome!" << std::endl;
					break;
				}

				if (!cur_rname.size()) {			// First refernece name
					cur_rname = rname;
					continue;
				}

				// Append collected references.
				logger << Time << "FastaFile: `" << cur_rname << "' Num refs : " << refs->size() << std::endl;
				chromosomes[cur_rname] = refs;
				refs = new Chromosome;
				cur_rname = rname;
			}

			// Append references.
			unsigned start = static_cast<unsigned>(refs->size());
			for (unsigned i = 0; i < line.size(); ++i) {
				unsigned char ch = (line[i] == 'A') ? 1 : (line[i] == 'C') ? 2 : (line[i] == 'G') ? 4 : (line[i] == 'T') ? 8 : 0;
				if ((start + i) % 2 == 0)
					refs->push_back(ch);
				else
					(*refs)[refs->size() - 1] |= ch << 4;
			}
		}

		// Append last collected chromosome.
		if (refs->size()) {
			logger << Time << "FastaFile: `" << cur_rname << "' Num refs : " << refs->size() << std::endl;
			chromosomes[cur_rname] = refs;
		}
		else
			delete refs;
		return true;
	}

	bool ReadBinaryFile(const std::string& fasta_path)
	{
		logger << Time << "FastaFile: Reading binary input file -> " << fasta_path << std::endl;
		std::ifstream fasta_file(fasta_path);
		if (!fasta_file.is_open()) {
			logger << Time << "FastaFile: Could not read binary file!" << std::endl;
			return false;
		}

		uint32_t sig;
		fasta_file.read(reinterpret_cast<char*>(&sig), sizeof(sig));
		if (sig != BINARY_FILE_SIGNATURE) {
			logger << Time << "FastaFile: Signature is not correct for binary file!" << std::endl;
			return false;
		}

		uint32_t num_chromosomes = 0;
		fasta_file.read(reinterpret_cast<char*>(&num_chromosomes), sizeof(num_chromosomes));
		for (unsigned i = 0; i < num_chromosomes; ++i) {
			std::string ref_name = ReadRefName(fasta_file);
			uint32_t num_bytes;
			fasta_file.read(reinterpret_cast<char*>(&num_bytes), sizeof(num_bytes));
			Chromosome* refs = new Chromosome;
			refs->resize(num_bytes);
			fasta_file.read(reinterpret_cast<char*>(&(*refs)[0]), num_bytes);
			chromosomes[ref_name] = refs;
		}
		return true;
	}

	std::string ReadRefName(std::ifstream& fasta_file) const
	{
		char ch;
		std::string ref_name;
		fasta_file.read(&ch, 1);
		while (ch) {
			ref_name += ch;
			fasta_file.read(&ch, 1);
		}
		return ref_name;
	}

	char GetRefChar(const Chromosome& vect, const std::string& ref_name, int pos) const
	{
		--pos;
		const int idx = pos / 2;
		if (idx > static_cast<int>(vect.size())) {
			logger << Time << "WARNING: Could not find `" << pos << "' at `" << ref_name
				<< "' chromosome (MAX:" << vect.size() << ')' << std::endl;
			return '\0';
		}

		const unsigned char ch = (pos % 2 == 0 ? vect[idx] : (vect[idx] >> 4)) & 0x0f;
		const char ref = ch == 1 ? 'A' : ch == 2 ? 'C' : ch == 4 ? 'G' : ch == 8 ? 'T' : 'N';
		return ref;
	}

	ChromeMap chromosomes;

	// Lazy fields
	std::string fasta_path;
	std::string cur_name;
	Chromosome lazy_chromosome;
};

#endif

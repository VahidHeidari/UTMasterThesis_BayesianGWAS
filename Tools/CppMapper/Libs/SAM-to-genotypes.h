#ifndef SAM_TO_GENOTYPES_H_
#define SAM_TO_GENOTYPES_H_

#include <climits>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#include "fasta-file.h"

#ifndef USE_LAZY_GET_NUCLEOTIDE
#define FASTA_CONST		const
#else
#define FASTA_CONST
#endif



enum NUCLEOTIDE_FLAGS {
	A = 0x01,
	C = 0x02,
	G = 0x04,
	T = 0x08,
	N = 0x10,
	MATCHED = 0x10000000
};



constexpr bool IS_BINARY_OUTPUT = true;
constexpr bool IS_COMPRESSED_OUTPUT = true;
constexpr std::ios_base::openmode OUT_FILE_FLAGS = (IS_BINARY_OUTPUT || IS_COMPRESSED_OUTPUT)
	? std::ios_base::binary | std::ios_base::out
	: std::ios_base::out;

constexpr int NUM_SAM_REC_PRINT = 500000;
constexpr int BITS_PER_GENOTYPE = 2;
constexpr int BUFFER_BITS = sizeof(uint32_t) * CHAR_BIT;



static inline uint32_t GetNucleotideFlag(char nucleotide)
{
	uint32_t flg = 0;
	switch (nucleotide) {
	case 'A': flg = NUCLEOTIDE_FLAGS::A; break;
	case 'C': flg = NUCLEOTIDE_FLAGS::C; break;
	case 'G': flg = NUCLEOTIDE_FLAGS::G; break;
	case 'T': flg = NUCLEOTIDE_FLAGS::T; break;
	case 'N': flg = NUCLEOTIDE_FLAGS::N; break;
	default: break;
	}
	return flg;
}

struct GenotypeRec
{
	inline bool HasNucleotide(int flag) const { return static_cast<int>(nucleotides & flag) == flag; }
	inline bool IsMatched() const { return HasNucleotide(NUCLEOTIDE_FLAGS::MATCHED); }

	inline bool HasNucleotide(char nucleotide) const
	{
		int flg = GetNucleotideFlag(nucleotide);
		return HasNucleotide(flg);
	}

	inline void AddNucleotide(char nucleotide)
	{
		const int flg = GetNucleotideFlag(nucleotide);
		nucleotides |= flg;
	}

	inline void SetCIGAR(char flag, const std::string& ref_name, FASTA_CONST FastaFile& fasta)
	{
		if (flag == 'X')
			return;

		if (flag == 'M') {
			char nuc = fasta.GetNucleotide(ref_name, pos);
			if (HasNucleotide(nuc))
				nucleotides |= NUCLEOTIDE_FLAGS::MATCHED;
		}
	}

	inline int GetNumFlags() const
	{
		return HasNucleotide(NUCLEOTIDE_FLAGS::A)
			+ HasNucleotide(NUCLEOTIDE_FLAGS::C)
			+ HasNucleotide(NUCLEOTIDE_FLAGS::G)
			+ HasNucleotide(NUCLEOTIDE_FLAGS::T);
	}

	inline bool IsBiallelic() const
	{
		const int num_flags = GetNumFlags();
		return num_flags == 0 || nucleotides == NUCLEOTIDE_FLAGS::N ? false : num_flags <= 2;
	}

	inline int GetGenotype() const
	{
		const int num_flags = GetNumFlags();
		if (num_flags == 2)
			return 1;

		if (num_flags == 1 && IsMatched())
			return 0;

		return 2;
	}

	uint32_t pos;
	uint32_t nucleotides;
};

struct PrintFunctor
{
	PrintFunctor(std::ostream& os = std::cout) : os(os) {}

	inline void WriteRec(const GenotypeRec& r, char op, std::integral_constant<bool, true>)
	{
		const unsigned char op_and_nucleotides = (op == ' ' ? 0 : 0x80) | (r.nucleotides & 0x7F);
		os.write(reinterpret_cast<const char*>(&r.pos), sizeof(r.pos));
		os.write(reinterpret_cast<const char*>(&op_and_nucleotides), 1);
	}

	inline void WriteRec(const GenotypeRec& r, char op, std::integral_constant<bool, false>)
	{
		os << op << r.pos << '\t';
		if (r.HasNucleotide(A))
			os << 'A';
		if (r.HasNucleotide(C))
			os << 'C';
		if (r.HasNucleotide(G))
			os << 'G';
		if (r.HasNucleotide(T))
			os << 'T';
		if (r.HasNucleotide(N))
			os << 'N';
		os << std::endl;
	}

	inline void operator()(const GenotypeRec& r, char op)
	{
		WriteRec(r, op, std::integral_constant<bool, IS_BINARY_OUTPUT>());
	}

	std::ostream& os;
};

class GenotypeList
{
public:
	GenotypeList()
	{
		genotypes.reserve(500000000);		// Allocate static memory.
	}

	int GetNumLoci() const { return static_cast<int>(genotypes.size() + inserts.size()); }
	int GetMemorySize() const { return GetNumLoci() * sizeof(GenotypeRec); }

	void Reset()
	{
		genotypes.clear();
		inserts.clear();
	}

	void Add(int sam_pos, char nucleotide, char CIGAR_flag, const std::string& ref_name, FASTA_CONST FastaFile& fasta)
	{
		int idx = GetGenotypeIdxBinSearch(sam_pos);
		if (idx == static_cast<int>(genotypes.size())) {
			GenotypeRec new_rec = { static_cast<uint32_t>(sam_pos), GetNucleotideFlag(nucleotide) };
			new_rec.SetCIGAR(CIGAR_flag, ref_name, fasta);
			genotypes.push_back(new_rec);
			return;
		}

		if (static_cast<int>(genotypes[idx].pos) == sam_pos) {
			genotypes[idx].AddNucleotide(nucleotide);
			genotypes[idx].SetCIGAR(CIGAR_flag, ref_name, fasta);
			return;
		}

		GenotypeRec new_rec = { static_cast<uint32_t>(sam_pos), GetNucleotideFlag(nucleotide) };
		new_rec.SetCIGAR(CIGAR_flag, ref_name, fasta);
		genotypes.insert(genotypes.begin() + idx, new_rec);
	}

	void AddInsert(int sam_pos, int insert_idx, char nucleotide)
	{
		--sam_pos;
		int idx = GetGenotypeIdx(sam_pos);
		if (idx == static_cast<int>(inserts.size())) {
			GenotypeRec new_rec = { static_cast<uint32_t>(sam_pos), GetNucleotideFlag(nucleotide) };
			inserts.push_back(new_rec);
			return;
		}

		// Move to insert_idx, extend if needed.
		for (int i = 0; i < insert_idx; ++i, ++idx) {
			if (static_cast<int>(inserts.size()) == idx || static_cast<int>(inserts[idx].pos) != sam_pos) {
				GenotypeRec new_rec = { static_cast<uint32_t>(sam_pos), 0 };
				inserts.insert(inserts.begin() + idx, new_rec);
			}
		}

		// Check last index, and extend if needed.
		if (static_cast<int>(inserts.size()) == idx || static_cast<int>(inserts[idx].pos) != sam_pos) {
			GenotypeRec new_rec = { static_cast<uint32_t>(sam_pos), 0 };
			inserts.insert(inserts.begin() + idx, new_rec);
		}

		// Append nucleotide at current index.
		inserts[idx].AddNucleotide(nucleotide);
	}

	template <class Functor = PrintFunctor>
	void PrintListsCombined(Functor functor = PrintFunctor()) const
	{
		if (!inserts.size()) {
			for (const auto& r : genotypes)
				functor(r, ' ');
			return;
		}

		if (!genotypes.size()) {
			for (const auto& r : inserts)
				functor(r, '*');
			return;
		}

		size_t i_idx = 0, g_idx = 0;
		while (true) {
			if (i_idx >= inserts.size()) {
				while (g_idx < genotypes.size())
					functor(genotypes[g_idx++], ' ');
				return;
			}

			if (g_idx >= genotypes.size()) {
				while (i_idx < inserts.size())
					functor(inserts[i_idx++], '*');
				return;
			}

			while (g_idx < genotypes.size() && genotypes[g_idx].pos <= inserts[i_idx].pos)
				functor(genotypes[g_idx++], ' ');
			if (g_idx >= genotypes.size())
				continue;

			while (i_idx < inserts.size() && inserts[i_idx].pos < genotypes[g_idx].pos)
				functor(inserts[i_idx++], '*');
		}
	}

	void DumpCompressed(std::ostream& out_file) const
	{
		// Dump genotypes.
		for (unsigned i = 0; i < genotypes.size(); ++i) {
			const GenotypeRec& rec = genotypes[i];
			if (!rec.IsBiallelic())
				continue;						// Process only bi-allelic loci.

			// Find continuous segment length.
			uint32_t len = 1;
			for (unsigned j = i + 1; j < genotypes.size(); ++j) {
				if (genotypes[j].pos != rec.pos + len || !genotypes[j].IsBiallelic())
					break;

				++len;
			}

			// Write to file.
			out_file.write(reinterpret_cast<const char*>(&rec.pos), sizeof(rec.pos));
			out_file.write(reinterpret_cast<const char*>(&len), sizeof(len));
			int bit_pos = 0;
			uint32_t buff = 0;
			int num_writes = 0;
			for (unsigned j = 0; j < len; ++j) {
				if (bit_pos == BUFFER_BITS) {
					// Flush buffer.
					out_file.write(reinterpret_cast<const char*>(&buff), sizeof(buff));
					bit_pos = buff = 0;
					++num_writes;
				}

				int geno = genotypes[i + j].GetGenotype() << bit_pos;
				buff |= geno;
				bit_pos += BITS_PER_GENOTYPE;
			}

			// Flush last bytes.
			if (bit_pos) {
				out_file.write(reinterpret_cast<const char*>(&buff), sizeof(buff));
				++num_writes;
			}

			if (num_writes != ((len - 1 + 16) / 16))
				logger << Time << warning << "------ num_writes:" << num_writes
					<< "   len:" << len << std::endl;

			i += len - 1;			// Advance index.
		}

		const uint32_t EOF_POS = 0;
		out_file.write(reinterpret_cast<const char*>(&EOF_POS), sizeof(EOF_POS));
		out_file.write(reinterpret_cast<const char*>(&EOF_POS), sizeof(EOF_POS));
	}

	void Dump(std::ostream& out_file) const
	{
		if (IS_BINARY_OUTPUT) {
			uint32_t sz = static_cast<uint32_t>(genotypes.size());
			out_file.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
			sz = static_cast<uint32_t>(inserts.size());
			out_file.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
		}

		if (IS_COMPRESSED_OUTPUT)
			DumpCompressed(out_file);
		else {
			PrintFunctor wf = PrintFunctor(out_file);
			PrintListsCombined(wf);
		}
		out_file.flush();
	}

private:
	int GetGenotypeIdxBinSearch(int sam_pos) const
	{
		if (!genotypes.size())
			return 0;

		if (genotypes.size() == 1)
			return genotypes[0].pos >= static_cast<uint32_t>(sam_pos) ? 0 : 1;

		int left = 0, right = static_cast<int>(genotypes.size()) - 1;
		while (left < right) {
			if (genotypes[right].pos == static_cast<uint32_t>(sam_pos))
				return right;

			if (genotypes[left].pos == static_cast<uint32_t>(sam_pos))
				return left;

			const int mid = left + (right - left) / 2;
			if (genotypes[mid].pos == static_cast<uint32_t>(sam_pos))
				return mid;

			if (genotypes[mid].pos < static_cast<uint32_t>(sam_pos))
				left = mid + 1;
			else
				right = mid - 1;
		}

		return right > 0 && genotypes[right].pos > static_cast<uint32_t>(sam_pos) ? right : right + 1;
	}

	int GetGenotypeIdx(int sam_pos) const
	{
		if (!inserts.size())
			return 0;

		for (int i = 0; i < static_cast<int>(inserts.size()); ++i)
			if (inserts[i].pos >= static_cast<uint32_t>(sam_pos))
				return i;

		return static_cast<int>(inserts.size());
	}

	std::vector<GenotypeRec> genotypes;
	std::vector<GenotypeRec> inserts;
};



std::string GetHumanSize(long long sz);

void WriteToGenotype(const std::string& sam_path, const std::string& out_path,
		FASTA_CONST FastaFile& fasta);

void WriteBAMToGenotype(const std::string& bam_path, const std::string& out_path,
		FASTA_CONST FastaFile& fasta);

#endif

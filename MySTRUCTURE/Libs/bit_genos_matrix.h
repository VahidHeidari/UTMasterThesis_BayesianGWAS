#ifndef BIT_GENOS_MATRIX_H_
#define BIT_GENOS_MATRIX_H_

#include <climits>
#include <cstring>
#include <fstream>
#include <string>

//#define USE_VECTOR_GENOS		1

#ifdef USE_VECTOR_GENOS
#include <vector>
#endif

class BitGenosMatrix
{
public:
	typedef unsigned int GenotypeMatrixType;

	static constexpr int BITS_PER_GENOTYPE = 2;
	static constexpr int GENOTYPE_MASK = (1 << BITS_PER_GENOTYPE) - 1;
	static constexpr int INDIVS_PER_WORD = sizeof(GenotypeMatrixType) * CHAR_BIT / BITS_PER_GENOTYPE;;

	static inline int GetWordsPerRow(int num_indivs) { return (num_indivs + INDIVS_PER_WORD) / INDIVS_PER_WORD; }

	BitGenosMatrix()
		: num_indivs(0), num_loci(0), num_clusters(0), num_words_per_row(0)
#ifndef USE_VECTOR_GENOS
		, genos(nullptr)
#endif
	{}

	BitGenosMatrix(int num_indivs, int num_loci, int num_clusters)
		: num_indivs(num_indivs)
		, num_loci(num_loci)
		, num_clusters(num_clusters)
		, num_words_per_row(GetWordsPerRow(num_indivs))
#ifndef USE_VECTOR_GENOS
		, genos(nullptr)
#endif
	{
		Init();
	}

	~BitGenosMatrix()
	{
#ifndef USE_VECTOR_GENOS
		delete[] genos;
		genos = nullptr;
#endif
	}

	void Init(int num_indivs, int num_loci, int num_clusters)
	{
#ifdef USE_VECTOR_GENOS
		genos.clear();
#endif

		this->num_indivs = num_indivs;
		this->num_loci = num_loci;
		this->num_clusters = num_clusters;
		num_words_per_row = GetWordsPerRow(num_indivs);
		Init();
	}

	inline int GetNumIndivs() const { return num_indivs; }
	inline int GetNumLoci() const { return num_loci; }
	inline int GetNumClusters() const { return num_clusters; }

	inline int GetGeno(int indiv, int locus) const
	{
#ifdef USE_VECTOR_GENOS
		return genos[indiv][locus];
#else
		const int LOCUS_IDX = GetLocusIdx(locus);
		const int INDIV_IDX = GetIndivIdx(indiv);
		const int BIT_NUM = GetBitNum(indiv);
		const int GENO = genos[LOCUS_IDX + INDIV_IDX] >> BIT_NUM;
		return GENO & GENOTYPE_MASK;
#endif
	}

	inline void SetGeno(int indiv, int locus, int val)
	{
#ifdef USE_VECTOR_GENOS
		genos[indiv][locus] = static_cast<unsigned char>(val);
#else
		const int LOCUS_IDX = GetLocusIdx(locus);
		const int INDIV_IDX = GetIndivIdx(indiv);
		const int BIT_NUM = GetBitNum(indiv);
		const int GENO_IDX = LOCUS_IDX + INDIV_IDX;
		genos[GENO_IDX] &= ~(GENOTYPE_MASK << BIT_NUM);
		genos[GENO_IDX] |= val << BIT_NUM;
#endif
	}

	inline bool DumpText(std::string path) const
	{
		std::ofstream genos_file(path);
		if (!genos_file.is_open())
			return false;

		// Dump configs.
		genos_file << "NUM_INDIVS: " << GetNumIndivs() << std::endl;
		genos_file << "NUM_LOCI: " << GetNumLoci() << std::endl;
		genos_file << "NUM_CLUSTERS: " << GetNumClusters() << std::endl;
		genos_file << std::endl;

		// Dump genotypes.
		for (int i = 0; i < GetNumIndivs(); ++i) {
			for (int l = 0; l < GetNumLoci(); ++l)
				genos_file << GetGeno(i, l);
			genos_file << std::endl;
		}
		genos_file.close();
		return true;
	}

	inline bool ReadFromTextFile(std::string path)
	{
#ifdef USE_VECTOR_GENOS
		genos.clear();
#endif

		std::ifstream genos_file(path);
		if (!genos_file.is_open())
			return false;

		// Read configs.
		int num_indivs, num_loci, num_clusters;
		std::string tmp_str;
		genos_file >> tmp_str >> num_indivs >> tmp_str >> num_loci >> tmp_str >> num_clusters;

		this->num_indivs = num_indivs;
		this->num_loci = num_loci;
		this->num_clusters = num_clusters;
		num_words_per_row = GetWordsPerRow(num_indivs);
		Init();

		// Read genotype.
		char geno;
		for (int i = 0; i < num_indivs; ++i)
			for (int l = 0; l < num_loci; ++l) {
				genos_file >> geno;
				if (geno == '\r' || geno == '\n')
					continue;

				const int G = (geno == '0' ? 0 : geno == '1' ? 1 : 2);
				SetGeno(i, l, G);
			}
		return true;
	}

private:
	inline int GetLocusIdx(int locus) const { return num_words_per_row * locus; }
	static inline int GetIndivIdx(int indiv) {	return indiv / INDIVS_PER_WORD; }
	inline int GetBitNum(int indiv) const { return (indiv % INDIVS_PER_WORD) * BITS_PER_GENOTYPE; }

	void Init()
	{
#ifdef USE_VECTOR_GENOS
		genos.resize(GetNumIndivs());
		for (int i = 0; i < GetNumIndivs(); ++i)
			genos[i].resize(GetNumLoci(), 0);
#else
		const int NUM_ITEMS = num_words_per_row * num_loci;
		genos = new GenotypeMatrixType[NUM_ITEMS];
		memset(genos, 0, NUM_ITEMS * sizeof(GenotypeMatrixType));
#endif
	}

	int num_indivs;
	int num_loci;
	int num_clusters;
	int num_words_per_row;

#ifdef USE_VECTOR_GENOS
	std::vector<std::vector<unsigned char>> genos;
#else
	GenotypeMatrixType* genos;
#endif
};
#endif

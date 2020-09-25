#ifndef PARAMS_H_
#define PARAMS_H_

#include <fstream>
#include <vector>

#define BASE_PATH					"E:\\C++\\MySTRUCTURE\\"
#define DEF_NO_ADMIX_CONF			BASE_PATH "no_admix_conf.txt"
#define DEF_ADMIX_CONF				BASE_PATH "admix_conf.txt"
#define DEF_EXH_MOTAHARI			BASE_PATH "exh_motahari_conf.txt"
#define UPDATE_FREQ					5
#define LAMBDA						1.0
#define ALPHA						1.0
#define INFECTION_PROB_THRESHOLD	0.7

/// Change value for logging full site probibilities.
#define LOG_FULL_SITE_PROBS			false
#define LOG_BOTH_SITE_PROBS			true

/// Uncomment for logging entropy calculations.
//#define TEST_ENTROPY		1

/// Uncomment for ignoring 1 probability sites.
#define IGNORE_SITES_WITH_1_PORBS	1

typedef std::vector<unsigned char> Chromosome;
typedef std::vector<Chromosome> Genotype;
typedef std::vector<Genotype> Individuals;

class Params
{
public:
	enum
	{
		CONTROL_LABEL = 0,
		CASE_LABEL = 1
	};

	Params();

	bool Init(const char* input_path) { return Init(input_path, false); }
	bool InitExhMotahari(const char* input_path) { return Init(input_path, true); }

	/// Parameters
	int GetNumIndividuals() const { return num_indivs; }
	int GetNumLoci() const { return num_loci; }
	int GetNumChromosomes() const { return num_chromosomes; }
	int GetNumClusters() const { return num_clusters; }
	int GetNumIterations() const { return num_iters; }
	int GetNumBurnins() const { return num_burnins; }

	const Individuals& GetIndividuals() const { return individuals; }
	const Genotype& GetGenotype(int i) const { return GetIndividuals()[i]; }
	const Chromosome& GetChromosome(int i, int c) const { return GetGenotype(i)[c]; }
	int GetAllele(int i, int c, int l) const { return GetChromosome(i, c)[l]; }
	int GetLabel(int i) const { return labels[i]; }
	bool IsCaseIndividual(int i) const { return GetLabel(i) == CASE_LABEL; }
	bool IsControlIndividual(int i) const { return GetLabel(i) == CONTROL_LABEL; }

	/// Allele informations
	int GetNumUniqueAlleles() const { return num_unique_alleles; }
	int GetNumAlleles(int l) const { return num_alleles[l]; }
	const std::vector<int>& GetNumAlleles() const { return num_alleles; }
	int GetMaxAllele(int l) const { return max_allele[l]; }
	int GetMinAllele(int l) const { return min_allele[l]; }

	void Print() const;

private:
	void SetAllele(int i, int c, int l, int a) { individuals[i][c][l] = a; }
	void SetMaxAllele(int l, int allele) { max_allele[l] = allele; }
	void SetMinAllele(int l, int allele) { min_allele[l] = allele; }
	void SetLabel(int i, int label) { labels[i] = label; }

	void CollectGenotypeInfo();
	void AdjustAlleles();

	bool Init(const char* input_path, bool is_labeled);
	void ReadFixedParams(std::ifstream& input_file);
	void ReadAlleles(std::ifstream& input_file, int i, int c);

	/// Config parameters
	int num_indivs;
	int num_chromosomes;
	int num_loci;
	int num_clusters;
	int num_iters;
	int num_burnins;
	Individuals individuals;

	/// Genotypes information
	int num_unique_alleles;
	std::vector<int> num_alleles;
	std::vector<int> max_allele;
	std::vector<int> min_allele;

	/// Exhaustive Motahari labels (1:case, 0:control)
	std::vector<unsigned char> labels;
};

#endif

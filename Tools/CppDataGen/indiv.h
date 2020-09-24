#ifndef INDIV_H_
#define INDIV_H_

#include <vector>

class Indiv
{
public:
	typedef std::vector<int> Chromosome;
	typedef std::vector<Chromosome> Genotype;
	typedef std::vector<int> OriginTrack;

	Indiv();
	Indiv(const Indiv& other);

	static Indiv CreateNewIndiv(const Indiv& father, const Indiv& mother, int cluster_num);
	static void CreateHaploid(Indiv& indiv, const Indiv& father, const Indiv& mother);
	static void CreateDiploid(Indiv& indiv, const Indiv& father, const Indiv& mother);
	static bool IsAdmixed(const Indiv& father, const Indiv& mother);

	void Init(int num_loci, int cluster_num, int new_phenotype = 0, int num_chromosomes = 2);
	int GetNumChromosomes() const;
	int GetNumLoci() const;
	int GetAllele(int chromosome, int locus) const;
	void SetAllele(int chromosome, int locus, int allele);
	void Print() const;

	bool IsAdmixedCluster(int k) const;
	bool IsAdmixedChromosome() const;

	bool is_admixed;
	int phenotype;
	OriginTrack cluster;
	Genotype genotype;
};

#endif

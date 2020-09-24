#include "indiv.h"

#include <random>

#include "logger.h"

Indiv::Indiv()
	 : is_admixed(false)
	 , phenotype(0)
	 , cluster(0)
{}

Indiv::Indiv(const Indiv& other)
	: is_admixed(other.is_admixed)
	, phenotype(other.phenotype)
	, genotype(other.genotype)
	, cluster(other.cluster)
{}

Indiv Indiv::CreateNewIndiv(const Indiv& father, const Indiv& mother, int cluster_num)
{
	Indiv indiv;
	const int NUM_CHROMS = father.GetNumChromosomes();
	indiv.genotype.resize(NUM_CHROMS);
	indiv.cluster.resize(NUM_CHROMS);

	if (NUM_CHROMS >= 3)
		logger << warning << "Chromosome number greater than 2 is not implemented yet!" << std::endl;

	if (NUM_CHROMS == 1)
		CreateHaploid(indiv, father, mother);
	else
		CreateDiploid(indiv, father, mother);

	return indiv;
}

void Indiv::CreateHaploid(Indiv& indiv, const Indiv& father, const Indiv& mother)
{
	std::random_device dev;
	std::uniform_int<int> unif(0, 1);
	int c = unif(dev);
	if (c == 0) {
		indiv.is_admixed = father.is_admixed;
		indiv.genotype[0] = father.genotype[0];
		indiv.cluster[0] = father.cluster[0];
	} else {
		indiv.is_admixed = mother.is_admixed;
		indiv.genotype[0] = mother.genotype[0];
		indiv.cluster[0] = mother.cluster[0];
	}
}

void Indiv::CreateDiploid(Indiv& indiv, const Indiv& father, const Indiv& mother)
{
	std::random_device dev;
	std::uniform_int<int> unif(0, 1);
	int c = unif(dev);
	indiv.genotype[0] = father.genotype[c];
	indiv.cluster[0] = father.cluster[c];

	c = unif(dev);
	indiv.genotype[1] = mother.genotype[c];
	indiv.cluster[1] = mother.cluster[c];

	if (IsAdmixed(father, mother))
		indiv.is_admixed = true;
}

bool Indiv::IsAdmixed(const Indiv& father, const Indiv& mother)
{
	if (father.is_admixed || mother.is_admixed)
		return true;

	const bool is_chrom_admix = father.cluster != mother.cluster;
	return is_chrom_admix;
}

void Indiv::Init(int num_loci, int cluster_num, int new_phenotype, int num_chromosomes)
{
	phenotype = new_phenotype;
	genotype.clear();
	cluster.clear();
	genotype.resize(num_chromosomes);
	cluster.resize(num_chromosomes);
	for (int c = 0; c < num_chromosomes; ++c) {
		genotype[c].resize(num_loci);
		cluster[c] = cluster_num;
	}
}

int Indiv::GetNumChromosomes() const
{
	return genotype.size();
}

int Indiv::GetNumLoci() const
{
	return genotype.size() ? genotype[0].size() : 0;
}

int Indiv::GetAllele(int chromosome, int locus) const
{
	return genotype[chromosome][locus];
}

void Indiv::SetAllele(int chromosome, int locus, int allele)
{
	genotype[chromosome][locus] = allele;
}

void Indiv::Print() const
{
	for (int c = 0; c < GetNumChromosomes(); ++c) {
		for (int l = 0; l < GetNumLoci(); ++l) {
			logger << GetAllele(c, l);
			if (l + 1 < GetNumLoci())
				logger << "    ";
		}
		logger << std::endl;
	}
}

bool Indiv::IsAdmixedCluster(int k) const
{
	for (int c = 0; c < GetNumChromosomes(); ++c)
		if (cluster[c] != k)
			return true;

	return false;
}

bool Indiv::IsAdmixedChromosome() const
{
	const int first_cluster = cluster[0];
	for (int c = 1; c < GetNumChromosomes(); ++c)
		if (cluster[c] != first_cluster)
			return true;

	return false;
}

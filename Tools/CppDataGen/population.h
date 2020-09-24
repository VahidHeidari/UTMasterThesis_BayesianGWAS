#ifndef POPULATION_H_
#define POPULATION_H_

#include <vector>

#include "indiv.h"
#include "roulette-wheel.h"

class Frequencies;
class SimPars;

class Population
{
public:
	typedef std::vector<Indiv> Cluster;
	typedef std::vector<Cluster> Pops;

	typedef std::vector<RouletteWheel> RouletteWheelVect;
	typedef std::vector<RouletteWheelVect> AlleleRouletteWheel;

	void Init(const SimPars& sim_pars);
	void InitAdmix(const SimPars& sim_pars, const Population& population);

	int GetNumIndividuals(int cluster) const;
	int GetNumIndividuals() const;
	int GetNumClusters() const;
	int GetNumLoci() const;
	int GetNumChromosomes() const;
	const Indiv& GetIndividual(int cluster, int idx) const;
	const Indiv& GetRandomIndividual(int cluster) const;
	AlleleRouletteWheel CreateAlleleRouletteWheel(const Frequencies& freqs) const;
	void CreateNewGenotype(Indiv& indiv, AlleleRouletteWheel& arw) const;
	int GetAlleleCount(int cluster_num, int chromosome, int locus, int allele) const;
	Frequencies GetClusterFrequencies(int cluster_num, const SimPars& sim_pars) const;
	void PrintClustersFrequencies(const SimPars& sim_pars) const;
	void PrintIndividuals() const;
	void PrintAdmixProportions(const char* lable) const;
	void PrintNumAdmixIndivs() const;
	void StepSimulation();
	void InfectDisease(const SimPars& sim_pars);
	
	void CreateSTRUCTUREInput() const;
	void CreateMySTRUCTUREInput(const char* out_file_path) const;
	void CreateMySTRUCTUREInput() const;
	void CreateMySTRUCTUREAdmixInput() const;
	void CreateExhMotahariInput(const SimPars& sim_pars);
	void CreateBEAM3Input() const;

	Pops pops;
	std::vector<std::vector<int>> disease_labels;
};



class Parents
{
public:
	Parents();

	void SetParents(int father, int mother);
	Indiv BreedNewIndividual(const Population& pop, int cluster_num) const;

	int father_idx;
	int mother_idx;
};

#endif

#ifndef SIM_PARS_H_
#define SIM_PARS_H_

#include <fstream>
#include <string>
#include <vector>

#define BASE_PATH						"E:\\C++\\DataGen\\"
#define DEF_CONFIG_PATH					BASE_PATH "simulation_parameters.txt"
#define DEF_STRUCTURE_PATH				BASE_PATH "structure_input.txt"
#define DEF_MY_STRUCTURE_PATH			BASE_PATH "no_admix_conf.txt"
#define DEF_MY_STRUCTURE_ADMIX_PATH		BASE_PATH "admix_conf.txt"
#define DEF_EXH_MOTAHARI_PATH			BASE_PATH "exh_motahari_conf.txt"
#define DEF_BEAM3_PATH					BASE_PATH "beam3.txt"
#define EPSILON							0.0001

class SimPars;

class Frequencies
{
public:
	typedef std::vector<double> Freqs;
	typedef std::vector<Freqs> AlleleFreqs;
	typedef std::vector<AlleleFreqs> FreqVect;

	void Init(int num_loci, int num_alleles, int num_chromosomes = 2);
	void Init(const SimPars& sim_pars);

	int GetNumChromosomes() const;
	int GetNumLoci() const;
	int GetNumAlleles() const;
	double GetFreq(int chromosome, int locus, int allele) const;
	void SetFreq(int chromosome, int locus, int allele, double freq);
	void LogFrequencies(const SimPars& sim_pars) const;
	int GetMaxFreqAllele(int chromosome, int locus) const;

	FreqVect freqs;
};



/// Admixture parameters
struct AdmixPars
{
	int num_indivs;
	int num_admix;
	std::vector<int> clusters;
	std::vector<double> proportion;
};



/// Simulation parameters
class SimPars
{
public:
	bool Init(std::string config_path);

	double GetFreq(int cluster, int chromosome, int locus, int allele) const { return allele_freqs[cluster].GetFreq(chromosome, locus, allele); }
	void SetFreq(int cluster, int chromosome, int locus, int allele, double freq) { allele_freqs[cluster].SetFreq(chromosome, locus, allele, freq); }
	int GetMaxFreqAllele(int cluster, int chromosome, int locus) const { return allele_freqs[cluster].GetMaxFreqAllele(chromosome, locus); }
	const std::vector<int>& GetDiseaseAlleles(int cluster) const { return disease_allele[cluster]; }
	int GetInfectedLoci(int cluster, int i_th_loci) const { return infected_loci[cluster][i_th_loci]; }
	int GetDiseaseAllele(int cluster, int i_th_loci) const { return disease_allele[cluster][i_th_loci]; }
	bool IsInfectedLocus(int cluster, int locus) const;
	bool IsInfectedLocusInClusters(int locus) const;
	bool IsDiffMAFLocus(int locus) const;

	void ReadAlleleFrequencies(std::ifstream& sim_pars);
	void LogInputParameters() const;
	bool IsSimModeAdmixture() const;
	bool IsSimModeNoAdmixture() const;
	bool IsSimModeCaseControl() const;
	void ReadAdmixtureParams(std::ifstream& sim_pars);
	void ReadNumIndividualsInEachCluster(std::ifstream& sim_pars);
	void ReadCaseControlParams(std::ifstream& sim_pars);
	void InitCaseControlFrequencies();
	void AddMAFDifference();
	void InitDiseaseAllele();

	/// Fixed number simulation parameters
	std::string sim_mode;
	int num_chromosomes;
	int num_loci;
	int num_alleles;
	int iters;
	int num_clusters;

	/// Variable number simulation parameters
	std::vector<int> num_indivs;

	/// Allele frequencies
	std::vector<Frequencies> allele_freqs;

	/// Admixture parameters
	/// Fixed number parameters
	int admix_iters;
	int admix_clusters;

	/// Variable number parameters
	std::vector<AdmixPars> admix_pars;

	/// Case-Control parameters
	double diff_MAF_percent;
	int avg_infected_loci;
	std::vector<double> case_percent;
	std::vector<int> diff_loci;

	/// Disease model
	std::vector<std::vector<int>> infected_loci;
	std::vector<std::vector<int>> disease_allele;			/// Random generated disease allele combinations.
};

#endif

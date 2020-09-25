#ifndef ALLELE_FREQUENCIES_H_
#define ALLELE_FREQUENCIES_H_

#include <vector>

#include "dists.h"

typedef std::vector<double> AlleleFreqs;
typedef std::vector<AlleleFreqs> LocusFreqs;
typedef std::vector<LocusFreqs> Freqs;

class AlleleFrequencies
{
public:
	void Init(int clusters, int loci, const std::vector<int>& num_alleles);

	double GetAlleleFreq(int cluster, int locus, int allele) const { return freqs[cluster][locus][allele]; }
	void SetAlleleFreq(int cluster, int locus, int allele, double value) { freqs[cluster][locus][allele] = value; }

	const AlleleFreqs& GetAlleleFreqs(int cluster, int locus) const { return freqs[cluster][locus]; }
	AlleleFreqs& GetAlleleFreqs(int cluster, int locus) { return freqs[cluster][locus]; }

private:
	Freqs freqs;
};



typedef std::vector<int> AllelesCount;
typedef std::vector<AllelesCount> CountsOfLocus;
typedef std::vector<CountsOfLocus> AlleleCounts;

class AllelesCounts
{
public:
	void Init(int clusters, int loci, const std::vector<int>& num_alleles);

	void Increment(int cluster, int loci, int allele) { ++allele_counts[cluster][loci][allele]; }

	void ZeroAll();

	int GetNumClusters() const { return static_cast<int>(allele_counts.size()); }
	int GetNumLoci() const { return static_cast<int>(allele_counts[0].size()); }
	int GetNumAlleles(int k, int l) const { return static_cast<int>(allele_counts[k][l].size()); }
	int GetAlleleCount(int cluster, int locus, int allele) const { return allele_counts[cluster][locus][allele]; }

private:
	AlleleCounts allele_counts;
};



typedef std::vector<int> AlleleOrigin;
typedef std::vector<AlleleOrigin> ChromosomeOrigin;
typedef std::vector<ChromosomeOrigin> PopOrigin;

class AdmixZ
{
public:
	void Init(int num_indivs, int chromosomes, int loci, int clusters);

	const int GetOrigin(int i, int c, int l) const { return Z[i][c][l]; }
	void SetOrigin(int i, int c, int l, int new_k) { Z[i][c][l] = new_k; }

private:
	PopOrigin Z;
};



typedef std::vector<double> OriginProportion;
typedef std::vector<OriginProportion> IndividualsOrigin;

class AdmixProportions
{
public:
	void Init(int num_indivs, int clusters);

	const OriginProportion& GetOriginProportions(int i) const { return Q[i]; }
	OriginProportion& GetOriginProportions(int i) { return Q[i]; }
	double GetOriginProportion(int i, int k) const { return Q[i][k]; }
	void SetOriginProportion(int i, int k, int value) { Q[i][k] = value; }

private:
	IndividualsOrigin Q;
};



struct DiseaseSite
{
	int loci_start;
	int i_th_allele_comb;
	int i_th_loci_comb;
	int omega;
	int w;
	int num_cases;
	int num_contorls;
	double prob;
};

typedef std::vector<DiseaseSite> DiseaseSites;
typedef std::vector<DiseaseSites> Sites;

typedef std::vector<int> AllelesComb;
typedef std::vector<AllelesComb> Combs;

class DiseaseModel
{
public:
	void Init(int num_indivs, int num_loci, int num_clusters,
		int num_unique_alleles, int num_causal_sites, int epistasis_radius);

	void Update(int k, int loci_start, int i_th_allele_comb,
		int i_th_loci_comb, int omega, int w, int num_cases, int num_controls);

	void ClearSites();

	const AllelesComb& GetAlleleCombination(int i_th) const { return C[i_th]; }
	const Permuts& GetSiteCombination(int i_th) const { return sites_combs[i_th]; }
	int GetNumAllelesCombinations() const { return static_cast<int>(C.size()); }
	int GetNumSitesCombinations() const { return static_cast<int>(sites_combs.size()); }
	int GetEpistasisRadius() const { return epistasis_radius; }
	int GetNumLoci() const { return num_loci; }

	const DiseaseSite& GetDiseaseSite(int cluster, int loci) const { return sites[cluster][loci]; }

private:
	int GetNumCausalSites() const { return num_causal_sites; }
	int GetNumUniqueAlleles() const { return num_unique_alleles; }
	int GetNumClusters() const { return num_clusters; }
	DiseaseSite& GetDiseaseSite(int cluster, int locus) { return sites[cluster][locus]; }

	void GetCombination(int i_th, std::vector<int>& output) const;
	void InitSites();
	void InitAlleleCombinations();

	bool IsSiteProbabilityLess(int cluster, int start_locus, int i_th_loci_comb, double prob) const;

	void SetDiseaseSite(int cluster, int loci_start, int i_th_allele_comb,
		int i_th_loci_comb, int omega, int w, int num_cases,
		int num_controls, double prob);

	bool IsCombinationFileExist() const;
	void MakeCombinationFromFile();
	void DumpCombination() const;

	int num_causal_sites;		/// Average disease sites
	int num_unique_alleles;
	int epistasis_radius;
	int num_clusters;
	int num_loci;
	Sites sites;
	Combs C;
	PermutsVect sites_combs;
};

#endif

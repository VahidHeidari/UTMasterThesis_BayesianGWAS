#include "allele-frequencies.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <random>

#include "logger.h"
#include "params.h"

void AlleleFrequencies::Init(int clusters, int loci, const std::vector<int>& num_alleles)
{
	freqs.clear();
	freqs.resize(clusters);
	for (int c = 0; c < clusters; ++c) {
		freqs[c].resize(loci);
		for (int l = 0; l < loci; ++l)
			freqs[c][l].resize(num_alleles[l]);
	}
}



void AllelesCounts::Init(int clusters, int loci, const std::vector<int>& num_alleles)
{
	allele_counts.clear();
	allele_counts.resize(clusters);
	for (int k = 0; k < clusters; ++k) {
		allele_counts[k].resize(loci);
		for (int l = 0; l < loci; ++l)
			allele_counts[k][l].resize(num_alleles[l]);
	}
}

void AllelesCounts::ZeroAll()
{
	// Clear counters.
	for (int k = 0; k < GetNumClusters(); ++k)
		for (int l = 0; l < GetNumLoci(); ++l)
			for (int a = 0; a < GetNumAlleles(k, l); ++a)
				allele_counts[k][l][a] = 0;
}



void AdmixZ::Init(int num_indivs, int chromosomes, int loci, int clusters)
{
	std::random_device dev;
	std::uniform_int<int> unif(0, clusters - 1);

	Z.clear();
	Z.resize(num_indivs);
	for (int i = 0; i < num_indivs; ++i) {
		Z[i].resize(chromosomes);
		for (int c = 0; c < chromosomes; ++c) {
			Z[i][c].resize(loci);
			for (int l = 0; l < loci; ++l) {
				SetOrigin(i, c, l, unif(dev));
			}
		}
	}
}



void AdmixProportions::Init(int num_indivs, int clusters)
{
	Q.clear();
	Q.resize(num_indivs);
	for (int i = 0; i < num_indivs; ++i) {
		Q[i].resize(clusters);
		for (int c = 0; c < clusters; ++c)
			Q[i][c] = 0;
	}
}



void DiseaseModel::Init(int num_indivs, int num_loci, int num_clusters,
	int num_unique_alleles, int num_causal_sites, int epistasis_radius)
{
	this->num_causal_sites = num_causal_sites;
	this->num_unique_alleles = num_unique_alleles;
	this->epistasis_radius = std::min(epistasis_radius, num_loci);
	this->num_clusters = num_clusters;
	this->num_loci = num_loci;

	if (this->epistasis_radius != epistasis_radius)
		logger << warning << "Epistasis radius of " << epistasis_radius
				<< " could not be used! Instead " << GetEpistasisRadius()
				<< " is used because number of loci is " << num_loci << " at most." << std::endl;

	InitSites();

	sites_combs.clear();
	if (IsCombinationFileExist())
		MakeCombinationFromFile();
	else {
		FastMakeCombinations(GetEpistasisRadius(), GetNumCausalSites(), sites_combs);
		DumpCombination();
	}

	InitAlleleCombinations();
}

static bool IsNaN(double d)
{
	return d != d;
}

void DiseaseModel::Update(int k, int loci_start, int i_th_allele_comb,
	int i_th_loci_comb, int omega, int w, int num_cases, int num_controls)
{
	if (omega == 0)
		return;

	const double INFECTED_PROB = std::log((double)w / (double)omega);
	const double INFECTED_PROB_COMP = std::log(1.0 - INFECTED_PROB);
	const double LOG_PROB = num_cases * INFECTED_PROB + num_controls * INFECTED_PROB_COMP;
	const double PROB = std::exp(LOG_PROB);
	if (IsNaN(PROB) || PROB < INFECTION_PROB_THRESHOLD)
		return;

#ifdef IGNORE_SITES_WITH_1_PORBS
	if (PROB >= 1)
		return;
#endif

	if (IsSiteProbabilityLess(k, loci_start, i_th_loci_comb, PROB))
		SetDiseaseSite(k, loci_start, i_th_allele_comb, i_th_loci_comb, omega, w, num_cases, num_controls, PROB);
}

void DiseaseModel::ClearSites()
{
	for (int k = 0; k < GetNumClusters(); ++k)
		for (int l = 0; l < GetNumLoci(); ++l)
			memset(&sites[k][0], 0, sites[k].size() * sizeof(DiseaseSite));
}

void DiseaseModel::GetCombination(int i_th, std::vector<int>& output) const
{
	int cnt = i_th;
	for (int i = 0; i < GetNumCausalSites(); ++i) {
		const int idx = cnt % GetNumUniqueAlleles();
		output[i] = idx;
		cnt /= GetNumUniqueAlleles();
	}
}

void DiseaseModel::InitSites()
{
	sites.clear();
	sites.resize(GetNumClusters());
	for (int k = 0; k < GetNumClusters(); ++k) {
		for (int l = 0; l < GetNumLoci(); ++l) {
			DiseaseSite site;
			sites[k].push_back(site);
		}
	}
	ClearSites();
}

void DiseaseModel::InitAlleleCombinations()
{
	const int NUM_COMBS = (int)std::pow(GetNumUniqueAlleles(), GetNumCausalSites());
	C.clear();
	C.resize(NUM_COMBS);
	for (int i = 0; i < NUM_COMBS; ++i) {
		C[i].resize(GetNumCausalSites());
		GetCombination(i, C[i]);
	}
}

bool DiseaseModel::IsSiteProbabilityLess(int cluster, int start_locus, int i_th_loci_comb, double prob) const
{
	const Permuts& loci_comb = GetSiteCombination(i_th_loci_comb);
	for (int comb : loci_comb) {
		const DiseaseSite& site = GetDiseaseSite(cluster, comb + start_locus);
		if (site.prob > prob)
			return false;
	}

	return true;
}

void DiseaseModel::SetDiseaseSite(int cluster, int loci_start, int i_th_allele_comb,
	int i_th_loci_comb, int omega, int w, int num_cases, int num_controls, double prob)
{
	const Permuts& site_comb = GetSiteCombination(i_th_loci_comb);
	for (int comb : site_comb) {
		DiseaseSite& site = GetDiseaseSite(cluster, loci_start + comb);
		site.loci_start = loci_start;
		site.i_th_allele_comb = i_th_allele_comb;
		site.i_th_loci_comb = i_th_loci_comb;
		site.omega = omega;
		site.w = w;
		site.num_cases = num_cases;
		site.num_contorls = num_controls;
		site.prob = prob;
	}
}

static std::string GetCombionationFileName(int N, int R)
{
	std::ostringstream  ss_file_name;
	ss_file_name << BASE_PATH << "comb-" << N << '-' << R << ".txt";
	return ss_file_name.str();
}

bool DiseaseModel::IsCombinationFileExist() const
{
	std::string file_name = GetCombionationFileName(GetEpistasisRadius(), GetNumCausalSites());
	std::ifstream input_file(file_name);
	if (!input_file.is_open())
		return false;

	int n, r;
	input_file >> n >> r;
	if (n != GetEpistasisRadius() && r != GetNumCausalSites())
		return false;

	return true;
}

void DiseaseModel::MakeCombinationFromFile()
{
	std::string file_name = GetCombionationFileName(GetEpistasisRadius(), GetNumCausalSites());
	std::ifstream input_file(file_name);

	sites_combs.clear();

	int n, r, num_recs;
	input_file >> n >> r >> num_recs;
	for (int i = 0; i < num_recs; ++i) {
		Permuts perm(r);
		for (int j = 0; j < r; ++j) {
			int item;
			input_file >> item;
			perm[j] = item;
		}
		sites_combs.push_back(perm);
	}
}

void DiseaseModel::DumpCombination() const
{
	std::string file_name = GetCombionationFileName(GetEpistasisRadius(), GetNumCausalSites());
	std::ofstream output_file(file_name);
	output_file << GetEpistasisRadius() << ' ' << GetNumCausalSites() << std::endl;
	output_file << std::endl;
	output_file << sites_combs.size() << std::endl;
	for (const auto& perm : sites_combs) {
		for (int i : perm)
			output_file << i << ' ';
		output_file << std::endl;
	}
}

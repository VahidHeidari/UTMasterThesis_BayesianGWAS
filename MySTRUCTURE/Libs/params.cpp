#include "params.h"

#include <algorithm>
#include <limits>
#include <set>

#include "logger.h"

#define MAX_INDIVS	4
#define MAX_LOCI	6

Params::Params()
	: num_indivs(0)
	, num_chromosomes(0)
	, num_loci(0)
	, num_clusters(0)
	, num_iters(0)
	, num_burnins(0)
	, num_unique_alleles(0)
{}

void Params::Print() const
{
	logger << " Num Individuals    : " << GetNumIndividuals() << std::endl;
	logger << " Num Chromosomes    : " << GetNumChromosomes() << std::endl;
	logger << " Num Loci           : " << GetNumLoci() << std::endl;
	logger << " Num Clusters       : " << GetNumClusters() << std::endl;
	logger << " Num Iterations     : " << GetNumIterations() << std::endl;
	logger << " Num Burnins        : " << GetNumBurnins() << std::endl;
	logger << " Num Unique Alleles : " << GetNumUniqueAlleles() << std::endl;
	logger << std::endl;

	// Print some individual infos.
	logger << " Some Individuals : " << std::endl;
	for (int i = 0; i < std::min(MAX_INDIVS, GetNumIndividuals()); ++i) {
		const int i_idx = i < MAX_INDIVS / 2 ? i : (GetNumIndividuals() - MAX_INDIVS + i);
		logger << "  #" << (i_idx + 1) << std::endl;
		for (int c = 0; c < GetNumChromosomes(); ++c) {
			logger << "    ";
			for (int l = 0; l < std::min(MAX_LOCI, GetNumLoci()); ++l) {
				const int l_idx = l < MAX_LOCI / 2 ? l : (GetNumLoci() - MAX_LOCI + l);
				logger << GetAllele(i_idx, c, l_idx);
				if (l_idx + 1 == MAX_LOCI / 2)
					logger << "   ...   ";
				else if (l + 1 < GetNumLoci())
					logger << ' ';
			}
			logger << std::endl;
		}
		if (i + 1 == MAX_INDIVS / 2) {
			logger << "        ." << std::endl;
			logger << "        ." << std::endl;
			logger << "        ." << std::endl;
		}
	}
}

void Params::CollectGenotypeInfo()
{
	num_unique_alleles = 0;
	num_alleles.resize(GetNumLoci());
	max_allele.resize(GetNumLoci());
	min_allele.resize(GetNumLoci());
	std::set<int> all_alleles;

	for (int l = 0; l < GetNumLoci(); ++l) {
		int tmp_max_allele = std::numeric_limits<int>::min();
		int tmp_min_allele = std::numeric_limits<int>::max();
		num_alleles[l] = 0;
		std::set<int> locus_alleles;

		for (int i = 0; i < GetNumIndividuals(); ++i) {
			for (int c = 0; c < GetNumChromosomes(); ++c) {
				const int allele = GetAllele(i, c, l);

				// Check allele in current locus.
				const auto& loc_allele = locus_alleles.find(allele);
				if (loc_allele == locus_alleles.end()) {
					tmp_max_allele = std::max(tmp_max_allele, allele);
					tmp_min_allele = std::min(tmp_min_allele, allele);
					++num_alleles[l];
					locus_alleles.insert(allele);
				}

				// Check allele in all individuals.
				const auto& all_itr = all_alleles.find(allele);
				if (all_itr != all_alleles.end())
					continue;

				++num_unique_alleles;
				all_alleles.insert(allele);
			}
		}
		SetMaxAllele(l, tmp_max_allele);
		SetMinAllele(l, tmp_min_allele);
	}
}

static int GetIndexOf(int allele, const std::vector<int>& unique_alleles)
{
	for (unsigned i = 0; i < unique_alleles.size(); ++i)
		if (allele == unique_alleles[i])
			return static_cast<int>(i);

	logger << warning << "Invalid allele at GetIndexOf(allele:" << allele
			<< ", unique_alleles:" << unique_alleles.size() << ')' << std::endl;
	logger << "    ";
	for (const auto& i : unique_alleles)
		logger << i << ' ';
	logger << std::endl;
	return -1;		// We should never reach here!
}

void Params::AdjustAlleles()
{
	for (int l = 0; l < GetNumLoci(); ++l) {
		// Collect unique alleles in a set.
		std::set<int> set_alleles;
		for (int i = 0; i < GetNumIndividuals(); ++i) {
			for (int c = 0; c < GetNumChromosomes(); ++c) {
				const int allele = GetAllele(i, c, l);
				set_alleles.insert(allele);
			}
		}

		// Adjust alleles of current locus by index of sorted list.
		std::vector<int> locus_alleles;
		locus_alleles.insert(locus_alleles.begin(), set_alleles.begin(), set_alleles.end());
		for (int i = 0; i < GetNumIndividuals(); ++i) {
			for (int c = 0; c < GetNumChromosomes(); ++c) {
				const int allele = GetAllele(i, c, l);
				const int idx = GetIndexOf(allele, locus_alleles);
				SetAllele(i, c, l, idx);
			}
		}
	}
}

bool Params::Init(const char* input_path, bool is_labeled)
{
	std::ifstream input_file(input_path);
	if (!input_file.is_open()) {
		logger << "Could not open input file '" << input_path << "'!" << std::endl;
		return false;
	}

	ReadFixedParams(input_file);		// Read fixed parameters.

	// Read variable parameters.
	labels.clear();
	labels.resize(GetNumIndividuals(), 0);
	individuals.clear();
	individuals.resize(GetNumIndividuals());
	for (int i = 0; i < GetNumIndividuals(); ++i) {
		individuals[i].resize(GetNumChromosomes());
		for (int c = 0; c < GetNumChromosomes(); ++c) {
			if (is_labeled) {
				// Read label, for each chromosome and it is setted by the last chromosome.
				char lbl;
				input_file >> lbl;
				SetLabel(i, lbl - '0');
			}

			// Read alleles of chromosome.
			individuals[i][c].resize(GetNumLoci());
			individuals[i].resize(GetNumChromosomes());
			ReadAlleles(input_file, i, c);
		}
	}

	CollectGenotypeInfo();
	AdjustAlleles();
	return true;
}

void Params::ReadFixedParams(std::ifstream& input_file)
{
	input_file >> num_indivs;
	input_file >> num_chromosomes;
	input_file >> num_loci;
	input_file >> num_clusters;
	input_file >> num_iters;
	input_file >> num_burnins;
}

void Params::ReadAlleles(std::ifstream& input_file, int i, int c)
{
	for (int l = 0; l < GetNumLoci(); ++l) {
		char allele;
		input_file >> allele;
		SetAllele(i, c, l, allele - '0');
	}
}

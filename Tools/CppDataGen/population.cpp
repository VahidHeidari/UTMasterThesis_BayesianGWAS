#include "population.h"

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>

#include "logger.h"
#include "sim-pars.h"

void Population::Init(const SimPars& sim_pars)
{
	pops.clear();
	pops.resize(sim_pars.num_clusters);
	for (int k = 0; k < sim_pars.num_clusters; ++k) {
		const Frequencies& freqs = sim_pars.allele_freqs[k];			// Fetch allele frequency of cluster.
		AlleleRouletteWheel arw = CreateAlleleRouletteWheel(freqs);		// Create cluster frequencey roulette wheel.
		for (int i = 0; i < sim_pars.num_indivs[k]; ++i) {
			Indiv indiv;
			indiv.Init(sim_pars.num_loci, k, 0, sim_pars.num_chromosomes);
			CreateNewGenotype(indiv, arw);
			pops[k].push_back(indiv);
		}
	}
}

void Population::InitAdmix(const SimPars& sim_pars, const Population& population)
{
	pops.clear();
	pops.resize(sim_pars.admix_clusters);
	for (int k = 0; k < sim_pars.admix_clusters; ++k) {
		const AdmixPars& pars = sim_pars.admix_pars[k];
		RouletteWheel rw(pars.proportion);
		for (int i = 0; i < pars.num_indivs; ++i) {
			const int pop_idx = rw.GetNext();
			const int cluster = pars.clusters[pop_idx];
			const Indiv& indiv = population.GetRandomIndividual(cluster);
			pops[k].push_back(indiv);
		}
	}
}

int Population::GetNumIndividuals(int cluster) const
{
	return pops[cluster].size();
}

int Population::GetNumIndividuals() const
{
	int all_indivs = 0;
	for (int k = 0; k < GetNumClusters(); ++k)
		all_indivs += GetNumIndividuals(k);
	return all_indivs;
}

int Population::GetNumClusters() const
{
	return pops.size();
}

int Population::GetNumLoci() const
{
	const Indiv& indiv = GetIndividual(0, 0);
	return indiv.GetNumLoci();
}

int Population::GetNumChromosomes() const
{
	for (int k = 0; k < GetNumClusters(); ++k)
		if (GetNumIndividuals(k))
			return GetIndividual(k, 0).GetNumChromosomes();

	return 0;
}

const Indiv& Population::GetIndividual(int cluster, int idx) const
{
	return pops[cluster][idx];
}

const Indiv& Population::GetRandomIndividual(int cluster) const
{
	const Cluster& pop = pops[cluster];
	std::random_device dev;
	std::uniform_int<int> unif(0, pop.size() - 1);
	int idx = unif(dev);
	return pop[idx];
}

Population::AlleleRouletteWheel Population::CreateAlleleRouletteWheel(const Frequencies& freqs) const
{
	AlleleRouletteWheel arw;
	for (int c = 0; c < freqs.GetNumChromosomes(); ++c) {
		RouletteWheelVect roulette_vect;
		for (int l = 0; l < freqs.GetNumLoci(); ++l) {
			std::vector<double> chances;
			for (int a = 0; a < freqs.GetNumAlleles(); ++a) {
				double f = freqs.GetFreq(c, l, a);
				chances.push_back(f);
			}
			RouletteWheel rw(chances);
			roulette_vect.push_back(rw);
		}
		arw.push_back(roulette_vect);
	}
	return arw;
}

void Population::CreateNewGenotype(Indiv& indiv, AlleleRouletteWheel& arw) const
{
	const int NUM_CHROMS = indiv.GetNumChromosomes();
	const int NUM_LOCI = indiv.GetNumLoci();
	for (int c = 0; c < NUM_CHROMS; ++c) {
		RouletteWheelVect& rwv = arw[c];
		for (int l = 0; l < NUM_LOCI; ++l) {
			RouletteWheel& rw = rwv[l];
			int allele = rw.GetNext();
			indiv.genotype[c][l] = allele;
		}
	}
}

int Population::GetAlleleCount(int cluster_num, int chromosome, int locus, int allele) const
{
	int count = 0;
	const Cluster& cluster = pops[cluster_num];
	for (int i = 0; i < GetNumIndividuals(cluster_num); ++i) {
		const Indiv& indiv = cluster[i];
		const int indiv_allele = indiv.GetAllele(chromosome, locus);
		if (allele == indiv_allele)
			++count;
	}
	return count;
}

Frequencies Population::GetClusterFrequencies(int cluster_num, const SimPars& sim_pars) const
{
	Frequencies freqs;
	freqs.Init(sim_pars.num_loci, sim_pars.num_alleles, sim_pars.num_chromosomes);
	for (int i = 0; i < GetNumIndividuals(cluster_num); ++i) {
		for (int c = 0; c < sim_pars.num_chromosomes; ++c) {
			for (int l = 0; l < sim_pars.num_loci; ++l) {
				for (int a = 0; a < sim_pars.num_alleles; ++a) {
					int count = GetAlleleCount(cluster_num, c, l, a);
					double frq = (double)count / (double)GetNumIndividuals(cluster_num);
					freqs.SetFreq(c, l, a, frq);
				}
			}
		}
	}
	return freqs;
}

void Population::PrintClustersFrequencies(const SimPars& sim_pars) const
{
#if _DEBUG
	for (int k = 0; k < GetNumClusters(); ++k) {
		logger << "Pop #" << k + 1 << "    size:" << GetNumIndividuals(k) << std::endl;
		Frequencies freqs = GetClusterFrequencies(k, sim_pars);
		freqs.LogFrequencies(sim_pars);
		logger << std::endl;
	}
#endif
}

void Population::PrintIndividuals() const
{
	int ind_idx = 1;
	for (int k = 0; k < GetNumClusters(); ++k) {
		for (int i = 0; i < GetNumIndividuals(k); ++i, ++ind_idx) {
			const Indiv& indiv = GetIndividual(k, i);
			logger << "Indiv (i:" << ind_idx << ", Cluster:" << k + 1 << ", Orig:[";
			for (int c = 0; c < GetNumChromosomes(); ++c) {
				logger << indiv.cluster[c] + 1;
				if (c + 1 < GetNumChromosomes())
					logger << ", ";
			}
			logger << ']';
			if (indiv.is_admixed || indiv.IsAdmixedChromosome())
				logger << ", IsAdmix:" << std::boolalpha << indiv.is_admixed;
			else if (indiv.IsAdmixedCluster(k))
				logger << ", AdmixedCluster";
			logger << ")" << std::endl;
			indiv.Print();
			logger << std::endl;
		}
		logger << std::endl << std::endl;
	}
}

void Population::PrintAdmixProportions(const char* lable) const
{
	logger << "Admix Proportions ";
	if (lable && strlen(lable))
		logger << "of '" << lable << "' ";
	//const auto prec = std::cout.precision();
	logger << ":    " ;//<< std::setprecision(2);

	for (int k = 0; k < GetNumClusters(); ++k) {
		double num_admix = 0.0;
		for (int i = 0; i < GetNumIndividuals(k); ++i) {
			const Indiv& indiv = GetIndividual(k, i);
			if (indiv.is_admixed)
				num_admix += 1.0;
		}
		const double NUM_INDIVS = (double)GetNumIndividuals(k);
		const double percent = num_admix / NUM_INDIVS * 100.0;
		logger << k + 1 << ':' << percent << '%';
		if (k + 1 < GetNumClusters())
			logger << "    ";
	}
	//logger << std::setprecision(prec);
	logger << std::endl;
}

void Population::PrintNumAdmixIndivs() const
{
	logger << "Num Admixed Indivs :     ";
	for (int k = 0; k < GetNumClusters(); ++k) {
		int num_admix = 0;
		for (int i = 0; i < GetNumIndividuals(k); ++i) {
			const Indiv& indiv = GetIndividual(k, i);
			if (indiv.is_admixed || indiv.IsAdmixedCluster(k) || indiv.IsAdmixedChromosome())
				++num_admix;
		}
		logger << k + 1 << ':' << num_admix << '/' << GetNumIndividuals(k);
		if (k + 1 < GetNumClusters())
			logger << "    ";
	}
	logger << std::endl;
}

static int GetMother(int father_idx, std::uniform_int<int>& unif, std::random_device& dev)
{
	int mother_idx;
	do {
		mother_idx = unif(dev);
	} while (mother_idx == father_idx);
	return mother_idx;
}

void Population::StepSimulation()
{
	std::random_device dev;
	for (int k = 0; k < GetNumClusters(); ++k) {
		const int NUM_INDIVS = GetNumIndividuals(k);
		Cluster new_pop;
		Parents parents;
		const Cluster& cluster = pops[k];
		std::uniform_int<int> unif(0, NUM_INDIVS - 1);

		while (new_pop.size() != NUM_INDIVS / 2) {							// Breed new individuals.
			const int father = unif(dev);
			const int mother = GetMother(father, unif, dev);
			parents.SetParents(father, mother);
			Indiv new_indiv = parents.BreedNewIndividual(*this, k);
			new_pop.push_back(new_indiv);
		}

		while (new_pop.size() != NUM_INDIVS) {								// Transfer previous generation individuals to new population.
			int idx = unif(dev);
			const Indiv& indiv = GetIndividual(k, idx);
			new_pop.push_back(indiv);
		}
		pops[k] = new_pop;
	}
}

static void WriteIndividualGenotype(std::ofstream& str_out, int label, const Indiv& indiv)
{
	for (int c = 0; c < indiv.GetNumChromosomes(); ++c) {
		str_out << label << "    ";
		for (int l = 0; l < indiv.GetNumLoci(); ++l) {
			int allele = indiv.GetAllele(c, l);
			str_out << allele;
			if (l + 1 < indiv.GetNumLoci())
				str_out << ' ';
		}		// Allele
		str_out << std::endl;
	}		// Chromosome
}

void Population::InfectDisease(const SimPars& sim_pars)
{
	std::random_device dev;

	disease_labels.clear();
	disease_labels.resize(GetNumClusters());

	for (int k = 0; k < GetNumClusters(); ++k) {
		disease_labels[k].resize(sim_pars.num_indivs[k], 0);
		std::uniform_int<int> unif(0, sim_pars.num_indivs[k] - 1);
		const int NUM_CASE = (int)std::ceil(sim_pars.case_percent[k] * sim_pars.num_indivs[k]);
		int num_case = 0;
		do {
			const int INDIV = unif(dev);
			if (disease_labels[k][INDIV])
				continue;

			++num_case;
			disease_labels[k][INDIV] = 1;
			Indiv& indiv = pops[k][INDIV];
			for (int i = 0; i < sim_pars.avg_infected_loci; ++i) {
				const int LOCI = sim_pars.GetInfectedLoci(k, i);
				const int ALLELE = sim_pars.GetDiseaseAllele(k, i);
				for (int c = 0; c < GetNumChromosomes(); ++c)
					indiv.SetAllele(c, LOCI, ALLELE);
			}
		} while (num_case < NUM_CASE);
	}
}

static std::string MakeExhMotahariFileName(const SimPars& sim_pars)
{
	std::ostringstream ss;
	ss << BASE_PATH << "exh_motahari_conf-L";
	for (int k = 0; k < sim_pars.num_clusters; ++k) {
		ss << '(';
		for (int i = 0; i < sim_pars.avg_infected_loci; ++i) {
			ss << sim_pars.GetInfectedLoci(k, i) + 1;
			if (i + 1 < sim_pars.avg_infected_loci)
				ss << '-';
		}
		ss << ')';
	}

	ss << "-A";
	for (int k = 0; k < sim_pars.num_clusters; ++k) {
		ss << '(';
		for (int i = 0; i < sim_pars.avg_infected_loci; ++i) {
			ss << sim_pars.GetDiseaseAllele(k, i);
			if (i + 1 < sim_pars.avg_infected_loci)
				ss << '-';
		}
		ss << ')';
	}

	ss << "-D(";
	for (unsigned i = 0; i < sim_pars.diff_loci.size(); ++i) {
		ss << sim_pars.diff_loci[i] + 1;
		if (i + 1 < sim_pars.diff_loci.size())
			ss << '-';
	}
	ss << ").txt";
	return ss.str();
}

void Population::CreateSTRUCTUREInput() const
{
	std::ofstream str_out(DEF_STRUCTURE_PATH);

	// Print first row: marker names.
	for (int i = 0; i < GetNumLoci(); ++i) {
		str_out << i + 1;
		if (i + 1 < GetNumLoci())
			str_out << ' ';
	}
	str_out << std::endl;

	int indiv_label = 1;
	for (int k = 0; k < GetNumClusters(); ++k) {
		const Cluster& cluster = pops[k];
		for (unsigned i = 0; i < cluster.size(); ++i) {
			const Indiv& indiv = cluster[i];		// Write genotype data of current individual.
			WriteIndividualGenotype(str_out, indiv_label, indiv);
			++indiv_label;
		}		// Individuals
	}		// Clusters
}

void Population::CreateMySTRUCTUREInput(const char* out_file_path) const
{
	std::ofstream out_file(out_file_path);
	out_file << GetNumIndividuals() << std::endl;
	out_file << GetNumChromosomes() << std::endl;
	out_file << GetNumLoci() << std::endl;
	out_file << GetNumClusters() << std::endl;
	out_file << "20000" << std::endl;
	out_file << "10000" << std::endl;
	for (int k = 0; k < GetNumClusters(); ++k) {
		for (int i = 0; i < GetNumIndividuals(k); ++i) {
			const Indiv& indiv = GetIndividual(k, i);
			for (int c = 0; c < indiv.GetNumChromosomes(); ++c) {
				for (int l = 0; l < indiv.GetNumLoci(); ++l) {
					out_file << indiv.GetAllele(c, l);
					if (l + 1 < indiv.GetNumLoci())
						out_file << ' ';
				}
				out_file << std::endl;
			}
		}
	}
}

void Population::CreateMySTRUCTUREInput() const
{
	CreateMySTRUCTUREInput(DEF_MY_STRUCTURE_PATH);
}

void Population::CreateMySTRUCTUREAdmixInput() const
{
	CreateMySTRUCTUREInput(DEF_MY_STRUCTURE_ADMIX_PATH);
}

void Population::CreateExhMotahariInput(const SimPars& sim_pars)
{
	std::string out_file_name = MakeExhMotahariFileName(sim_pars);
	std::ofstream out_file(out_file_name);
	out_file << GetNumIndividuals() << std::endl;
	out_file << GetNumChromosomes() << std::endl;
	out_file << GetNumLoci() << std::endl;
	out_file << GetNumClusters() << std::endl;
	out_file << "20000" << std::endl;
	out_file << "10000" << std::endl << std::endl;
	for (int k = 0; k < GetNumClusters(); ++k) {
		for (int i = 0; i < GetNumIndividuals(k); ++i) {
			const Indiv& indiv = GetIndividual(k, i);
			const int LABEL = disease_labels[k][i];
			for (int c = 0; c < indiv.GetNumChromosomes(); ++c) {
				out_file << LABEL << "        ";
				for (int l = 0; l < indiv.GetNumLoci(); ++l) {
					out_file << indiv.GetAllele(c, l);						// Output current allele.

					if (l + 1 < indiv.GetNumLoci()) {
						if (sim_pars.IsInfectedLocus(k, l + 1))
							out_file << "   ";								// Output spaces before infected loci.
						else if (sim_pars.IsInfectedLocus(k, l))
							out_file << "    ";								// Output spaces after infected loci.
						else
							out_file << ' ';
					}
					out_file.flush();
				}
				out_file << std::endl;
			}		// Chromosome
			out_file << std::endl;
		}		// Individual
		out_file << std::endl << std::endl;
	}		// Cluster
}

static int GetGenotype(const Indiv& indiv, int locus)
{
	static const int GENOTYPES[] =
	{
		// 0  1  2  3
		   0, 1, 0, 1,		// 0
		   1, 2, 1, 2,		// 1
		   0, 1, 0, 1,		// 2
		   1, 2, 1, 2,		// 3
	};

	const int ALLELES[] = { indiv.GetAllele(0, locus), indiv.GetAllele(1, locus) };
	const int IDX = ALLELES[0] * 4 + ALLELES[1];
	const int GENO = GENOTYPES[IDX];
	return GENO;
}

void Population::CreateBEAM3Input() const
{
	std::ofstream out_file(DEF_BEAM3_PATH);
	out_file << "ID Chr Pos ";
	for (int k = 0; k < GetNumClusters(); ++k) {
		for (int i = 0; i < GetNumIndividuals(k); ++i) {
			out_file << disease_labels[k][i];

			if (k + 1 < GetNumClusters())
				out_file << ' ';
			else if (i + 1 < GetNumIndividuals(k))
				out_file << ' ';
		}
	}
	out_file << std::endl;

	int rs_num = 1;
	for (int l = 0; l < GetNumLoci(); ++l, ++rs_num) {
		out_file << "rs" << rs_num << " Chr1 " << 10000 + rs_num << ' ';
		for (int k = 0; k < GetNumClusters(); ++k) {
			for (int i = 0; i < GetNumIndividuals(k); ++i) {
				const Indiv& indiv = GetIndividual(k, i);
				const int GENO = GetGenotype(indiv, l);
				out_file << GENO;

				if (k + 1 < GetNumClusters())
					out_file << ' ';
				else if (i + 1 < GetNumIndividuals(k))
					out_file << ' ';
			}
		}

		if (l + 1 < GetNumLoci())
			out_file << std::endl;
	}
}



Parents::Parents()
{}

void Parents::SetParents(int father, int mother)
{
	father_idx = father;
	mother_idx = mother;
}

Indiv Parents::BreedNewIndividual(const Population& pop, int cluster_num) const
{
	const Indiv& father = pop.GetIndividual(cluster_num, father_idx);
	const Indiv& mother = pop.GetIndividual(cluster_num, mother_idx);
	Indiv indiv = Indiv::CreateNewIndiv(father, mother, cluster_num);
	return indiv;
}

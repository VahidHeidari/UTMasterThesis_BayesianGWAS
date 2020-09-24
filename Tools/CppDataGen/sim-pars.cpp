#include "sim-pars.h"

#include <algorithm>
#include <random>

#include "logger.h"

void Frequencies::Init(int num_loci, int num_alleles, int num_chromosomes)
{
	freqs.clear();
	freqs.resize(num_chromosomes);
	for (int c = 0; c < num_chromosomes; ++c) {
		AlleleFreqs& af = freqs[c];
		af.resize(num_loci);
		for (int l = 0; l < num_loci; ++l) {
			Freqs& f = af[l];
			f.resize(num_alleles);
		}
	}
}

void Frequencies::Init(const SimPars& sim_pars)
{
	Init(sim_pars.num_loci, sim_pars.num_alleles, sim_pars.num_chromosomes);
}

int Frequencies::GetNumChromosomes() const
{
	return freqs.size();
}

int Frequencies::GetNumLoci() const
{
	return freqs[0].size();
}

int Frequencies::GetNumAlleles() const
{
	return freqs[0][0].size();
}

double Frequencies::GetFreq(int chromosome, int locus, int allele) const
{
	const AlleleFreqs& af = freqs[chromosome];
	const Freqs& f = af[locus];
	double freq = f[allele];
	return freq;
}

void Frequencies::SetFreq(int chromosome, int locus, int allele, double freq)
{
	AlleleFreqs& af = freqs[chromosome];
	Freqs& f = af[locus];
	f[allele] = freq;
}

void Frequencies::LogFrequencies(const SimPars& sim_pars) const
{
	for (int c = 0; c < sim_pars.num_chromosomes; ++c) {
		logger << "    CHROMOSOME #" << c + 1 << std::endl;
		for (int a = 0; a < sim_pars.num_alleles; ++a) {
			logger << "        ";
			for (int l = 0; l < sim_pars.num_loci; ++l) {
				double f = GetFreq(c, l, a);
				logger << f << "    ";
			}
			logger << std::endl;
		}
		logger << std::endl;
	}
}

int Frequencies::GetMaxFreqAllele(int chromosome, int locus) const
{
	int idx = 0;
	double max_freq = std::numeric_limits<double>::min();
	for (int a = 0; a < GetNumAlleles(); ++a) {
		const double FREQ = GetFreq(chromosome, locus, a);
		if (max_freq < FREQ) {
			max_freq = FREQ;
			idx = a;
		}
	}
	return idx;
}



bool SimPars::Init(std::string config_path)
{
	std::ifstream sim_pars(config_path);
	if (!sim_pars.is_open()) {
		logger << "Could not open parameters file!" << std::endl;
		return false;
	}

	// Read fixed number parameters.
	sim_pars >> sim_mode;
	sim_pars >> num_chromosomes;
	sim_pars >> num_loci;
	sim_pars >> num_alleles;
	sim_pars >> iters;
	sim_pars >> num_clusters;

	// Read mode specific parameters.
	if (IsSimModeCaseControl()) {
		ReadCaseControlParams(sim_pars);
		InitCaseControlFrequencies();
		AddMAFDifference();
		InitDiseaseAllele();
	} else if (IsSimModeAdmixture() || IsSimModeNoAdmixture()) {
		ReadNumIndividualsInEachCluster(sim_pars);						// Read number of individuals in each of clusters.
		ReadAlleleFrequencies(sim_pars);								// Read allele chances of each cluster.
		if (IsSimModeAdmixture())
			ReadAdmixtureParams(sim_pars);								// Read admixture parameters if any.
	}
	return true;
}

bool SimPars::IsInfectedLocus(int cluster, int locus) const
{
	return std::find(infected_loci[cluster].begin(), infected_loci[cluster].end(), locus) != infected_loci[cluster].end();
}

bool SimPars::IsInfectedLocusInClusters(int locus) const
{
	for (int k = 0; k < num_clusters; ++k)
		if (IsInfectedLocus(k, locus))
			return true;

	return false;
}

bool SimPars::IsDiffMAFLocus(int locus) const
{
	return std::find(diff_loci.begin(), diff_loci.end(), locus) != diff_loci.end();
}

void SimPars::ReadAlleleFrequencies(std::ifstream& sim_pars)
{
	double chance = 0.0;
	for (int k = 0; k < num_clusters; ++k) {
		Frequencies freq;
		freq.Init(*this);
		for (int c = 0; c < num_chromosomes; ++c) {
			for (int a = 0; a < num_alleles; ++a) {
				for (int l = 0; l < num_loci; ++l) {
					sim_pars >> chance;
					freq.SetFreq(c, l, a, chance);
				}
			}
		}
		allele_freqs.push_back(freq);
	}
}

void SimPars::LogInputParameters() const
{
	logger << "Simulation Parameters:" << std::endl;
	logger << "    SIM MODE        : " << sim_mode << std::endl;
	logger << "    NUM CHROMOSOMES : " << num_chromosomes << std::endl;
	logger << "    NUM LOCI        : " << num_loci << std::endl;
	logger << "    NUM ALLELES     : " << num_alleles << std::endl;
	logger << "    ITERS           : " << iters << std::endl;
	logger << "    NUM CLUSTERS    : " << num_clusters << std::endl;
	logger << "    NUM INDIVS      : " << num_indivs.size() << std::endl;

	// Print number of indivduals in each cluster.
	for (int k = 0; k < num_clusters; ++k)
		logger << "        " << num_indivs[k] << std::endl;
	logger << std::endl;

#if _DEBUG
	// Print allele frequencies in each cluster.
	logger << "    ALLELE FREQUENCIES :" << std::endl;
	for (int k = 0; k < num_clusters; ++k) {
		logger << "    CLUSTER #" << k << std::endl;
		allele_freqs[k].LogFrequencies(*this);
		logger << std::endl << std::endl;
	}
	logger << std::endl;
#endif

	if (IsSimModeCaseControl()) {
		logger << "   DIFF MAF PERCENT : " << diff_MAF_percent << std::endl;
		logger << "   AVG DISEASE LOCI : " << avg_infected_loci << std::endl;

		logger << "   CASE PERCENT     : " << std::endl;
		for (int k = 0; k < num_clusters; ++k) {
			logger << "      #" << k + 1 << ":" << case_percent[k] << '%';
			if (k + 1 == num_clusters)
				logger << std::endl;
			else
				logger << "    ";
		}

		logger << std::endl;
		logger << "    INFECTED LOCI   : " << std::endl;
		for (int k = 0; k < num_clusters; ++k) {
			logger  << "      #" << k + 1 << ":    ";
			for (int i = 0; i < avg_infected_loci; ++i) {
				const int INF_LOCUS = infected_loci[k][i] + 1;
				logger << INF_LOCUS;
				if (i + 1 < avg_infected_loci)
					logger << "  ";
			}
			logger << std::endl;
		}

		logger << std::endl;
		logger << "    DIFF MAF LOCI   : ";
		for (unsigned l = 0; l < diff_loci.size(); ++l) {
			logger << diff_loci[l] + 1;
			if (l + 1 < diff_loci.size())
				logger << ' ';
		}
		logger << std::endl;

		logger << std::endl;
		logger << "    DISEASE ALLELE  :" << std::endl;
		for (int k = 0; k < num_clusters; ++k) {
			logger << "      #" << k + 1 << " :     ";
			for (int i = 0; i < avg_infected_loci; ++i) {
				logger << GetDiseaseAllele(k, i);
				if (i + 1 < avg_infected_loci)
					logger << "  ";
			}
			logger << std::endl;
		}
	} else if (IsSimModeAdmixture()) {
		// Print admixture parameters.
		logger << "    ADMIX PARS:" << std::endl;
		logger << "        ADMIX ITERS    : " << admix_iters << std::endl;
		logger << "        ADMIX CLUSTERS : " << admix_clusters << std::endl;
		logger << "            CLUSTERS:" << std::endl;
		for (int i = 0; i < admix_clusters; ++i) {
			const AdmixPars& pars = admix_pars[i];
			logger << "               NUM INDIVS : " << pars.num_indivs << std::endl;
			logger << "               NUM ADMIX  : " << pars.num_admix << std::endl;
			logger << "               ";
			for (int j = 0; j < pars.num_admix; ++j) {
				logger << pars.clusters[j];
				if (j + 1 < pars.num_admix)
					logger << "   ";
			}
			logger << std::endl;

			logger << "               ";
			for (int j = 0; j < pars.num_admix; ++j) {
				logger << pars.proportion[j];
				if (j + 1 < pars.num_admix)
					logger << "   ";
			}
			logger << std::endl;

			if (i + 1 < admix_clusters)
				logger << std::endl;
		}
	}
}

bool SimPars::IsSimModeAdmixture() const
{
	return sim_mode == "ADMIXTURE";
}

bool SimPars::IsSimModeNoAdmixture() const
{
	return sim_mode == "NO_ADMIXTURE";
}

bool SimPars::IsSimModeCaseControl() const
{
	return sim_mode == "CASE_CONTROL";
}

void SimPars::ReadAdmixtureParams(std::ifstream& sim_pars)
{
	sim_pars >> admix_iters;
	sim_pars >> admix_clusters;

	for (int i = 0; i < admix_clusters; ++i) {
		AdmixPars pars;
		sim_pars >> pars.num_indivs;
		sim_pars >> pars.num_admix;

		int cluster = 0;
		for (int j = 0; j < pars.num_admix; ++j) {
			sim_pars >> cluster;
			pars.clusters.push_back(cluster);
		}

		double prop = 0.0;
		for (int j = 0; j < pars.num_admix; ++j) {
			sim_pars >> prop;
			pars.proportion.push_back(prop);
		}
		admix_pars.push_back(pars);
	}
}

void SimPars::ReadNumIndividualsInEachCluster(std::ifstream& sim_pars)
{
	int num_indiv = 0;
	for (int i = 0; i < num_clusters; ++i) {
		sim_pars >> num_indiv;
		num_indivs.push_back(num_indiv);
	}
}

void SimPars::ReadCaseControlParams(std::ifstream& sim_pars)
{
	sim_pars >> diff_MAF_percent;
	sim_pars >> avg_infected_loci;
	ReadNumIndividualsInEachCluster(sim_pars);

	// Read percent of case individuals in each cluster.
	case_percent.clear();
	for (int k = 0; k < num_clusters; ++k) {
		double per;
		sim_pars >> per;
		case_percent.push_back(per);
	}

	// Read infected loci position in each cluster.
	infected_loci.clear();
	infected_loci.resize(num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		for (int l = 0; l < avg_infected_loci; ++l) {
			int inf_locus;
			sim_pars >> inf_locus;
			infected_loci[k].push_back(inf_locus - 1);
		}
	}
}

static void FillProps(std::vector<double>& output, double sum = 1.0)
{
	std::random_device dev;
	std::uniform_real<double> unif(0, sum);

	std::vector<double> rnd_points(output.size() - 1);			// Select some random points between [0, 1].
	for (unsigned i = 0; i < rnd_points.size(); ++i)
		rnd_points[i] = unif(dev);
	rnd_points.push_back(sum);
	std::sort(rnd_points.begin(), rnd_points.end());

	double prev_point = 0.0;
	for (unsigned i = 0; i < output.size(); ++i) {				// Fill output frequencies.
		output[i] = rnd_points[i] - prev_point;
		prev_point = rnd_points[i];
	}
}

void SimPars::InitCaseControlFrequencies()
{
	std::vector<double> props(num_alleles, 0);

	allele_freqs.clear();

	Frequencies freqs;
	freqs.Init(*this);
	for (int l = 0; l < num_loci; ++l) {
		FillProps(props);
		for (int c = 0; c < num_chromosomes; ++c)
			for (int a = 0; a < num_alleles; ++a)
			freqs.SetFreq(c, l, a, props[a]);
	}
	allele_freqs.push_back(freqs);

	// Copy to all clusters.
	for (int k = 1; k < num_clusters; ++k)
		allele_freqs.push_back(freqs);
}

void SimPars::AddMAFDifference()
{
	const double LOCI_PROP = std::ceil(diff_MAF_percent * (double)num_loci);
	const unsigned NUM_MAF_LOCI = static_cast<int>(LOCI_PROP);
	std::random_device dev;
	std::uniform_int<int> unif(0, num_loci - 1);

	// Select random MAF different loci.
	diff_loci.clear();
	do {
		const int LOCUS = unif(dev);
		if (IsDiffMAFLocus(LOCUS) || IsInfectedLocusInClusters(LOCUS))
			continue;

		diff_loci.push_back(LOCUS);
	} while (diff_loci.size() < NUM_MAF_LOCI);
	std::sort(diff_loci.begin(), diff_loci.end());

	// Change frequencies.
	for (int l : diff_loci) {
		const int MAX_FREQ_ALLELE = GetMaxFreqAllele(0, 0, l);
		const double MAX_FREQ = allele_freqs[0].GetFreq(0, l, MAX_FREQ_ALLELE);
		for (int k = 0; k < num_clusters; ++k) {
			std::vector<double> props(num_alleles - 1);
			FillProps(props, 1.0 - MAX_FREQ);
			props.insert(props.begin() + MAX_FREQ_ALLELE, MAX_FREQ);
			for (int a = 0; a < num_alleles; ++a)
				for (int c = 0; c < num_chromosomes; ++c)
					SetFreq(k, c, l, a, props[a]);
		}
	}
}

void SimPars::InitDiseaseAllele()
{
	std::random_device dev;
	std::uniform_int<int> unif(0, num_alleles - 1);

	disease_allele.clear();
	disease_allele.resize(num_clusters);
	for (int k = 0; k < num_clusters; ++k)
		for (int i = 0; i < avg_infected_loci; ++i) {
			const int LOCUS = GetInfectedLoci(k, i);
			const int MAX_FREQ_ALLELE = GetMaxFreqAllele(k, 0, LOCUS);

			int inf_allele;
			do {
				inf_allele = unif(dev);
			} while (MAX_FREQ_ALLELE == inf_allele);
			disease_allele[k].push_back(inf_allele);
		}
}

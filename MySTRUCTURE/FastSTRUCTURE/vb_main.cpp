#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <random>
#include <vector>
#include <utility>

#include "bit_genos_matrix.h"
#include "logger.h"

// Comment to use raw pointers.
#define USE_VECTOR							1

// Uncomment only one of the following lines.
//#define READ_GENOTYPES_FROM_BINARY_FILE	1
//#define MAKE_RANDOM_FREQS					1

// Uncomment to calculate LLBO for checking convergence.
//#define USE_LLBO							1



typedef float FloatType;
typedef std::vector<std::vector<std::pair<double, double>>> FreqsVector;



static constexpr int NUM_INDIVS = 500;
static constexpr int NUM_LOCI = 100;
static constexpr int NUM_CLUSTERS = 2;

static constexpr int MAX_ITERS   = 10;
static constexpr int ITER_REPORT = 10;
static constexpr int MIN_ITERS = 5;

//static constexpr int LLBO_UPDATE = 10;
static constexpr double LLBO_EPSILON = 1;

static constexpr int NUM_CHROMOSOMES = 2;



static const std::string DUMP_PATH =
#if defined __linux__ || defined __CYGWIN__
	"";
#else
	"F:\\C++\\MySTRUCTURE\\";
#endif

static const std::string FREQS_PATH = DUMP_PATH + "freqs.txt";
static const std::string GENOS_PATH = DUMP_PATH + "genos.txt";



static double digammal(double x)
{
	double result = 0, xx, xx2, xx4;
	if (x <= 0) {
		logger << Time << " Input error!    digamma(x:" << x << " <= 0)" << std::endl;
		exit(1);
	}
	for (; x < 7; ++x)
		result -= 1 / x;
	x -= 1.0 / 2.0;
	xx = 1.0 / x;
	xx2 = xx*xx;
	xx4 = xx2*xx2;
	result += log(x) + (1. / 24.)*xx2 - (7.0 / 960.0)*xx4 + (31.0 / 8064.0)*xx4*xx2 - (127.0 / 30720.0)*xx4*xx4;
	return result;
}

inline static FloatType digamma(FloatType x)
{
	return static_cast<FloatType>(digammal(static_cast<double>(x)));
}

inline static double LogBeta(double x, double y)
{
	// beta function : (tgamma(x) * tgamma(y)) / tgamma(x + y)
	const double res = lgamma(x) + lgamma(y) - lgamma(x + y);
	return res;  
}

inline static double LogFact(int x)
{
	double res = 0.0;
	for (int i = 2; i < x; ++i)
		res += log(i);
	return res;
}

inline static double Comb(int n, int k)
{
	if (k == 0 || k == n)
		return 1.0;

	if (k == 1 || k == n - 1)
		return n;

	const double log_comb = LogFact(n) - LogFact(n - k) - LogFact(k);
	const double exp_comb = exp(log_comb);
	return exp_comb;
}



struct Z
{
	Z() : num_indivs(0), num_loci(0), num_clusters(0)
#ifndef USE_VECTOR
		, assignments(nullptr)
#endif
	{}

	~Z()
	{
#ifndef USE_VECTOR
		delete[] assignments; assignments = nullptr;
#endif
		num_indivs = num_loci = num_clusters = 0;
	}

	inline void Init(int num_indivs, int num_loci, int num_clusters)
	{
		this->num_indivs = num_indivs;
		this->num_loci = num_loci;
		this->num_clusters = num_clusters;

#ifdef USE_VECTOR
		assignments.resize(GetNumIndivs());
		for (int n = 0; n < GetNumIndivs(); ++n) {
			assignments[n].resize(NUM_CHROMOSOMES);
			for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
				assignments[n][c].resize(GetNumLoci());
				for (int l = 0; l < GetNumLoci(); ++l)
					assignments[n][c][l].resize(GetNumClusters());
			}
		}
#else
		const int NUM_ELEMS = GetNumIndivs() * GetNumLoci() * GetNumClusters() * NUM_CHROMOSOMES;
		assignments = new FloatType[NUM_ELEMS];
#endif
		// Initialize uniform.
		std::random_device rd;
		std::uniform_real_distribution<FloatType> uf(0, 1);
		for (int n = 0; n < GetNumIndivs(); ++n)
			for (int l = 0; l < GetNumLoci(); ++l)
				for (int k = 0; k < GetNumClusters(); ++k) {
					SetAssignment(n, 0, l, k, static_cast<FloatType>(1.0 + 0.1 * (0.5 - uf(rd))));
					SetAssignment(n, 1, l, k, static_cast<FloatType>(1.0 + 0.1 * (0.5 - uf(rd))));
				}
		Normalize();
	}

	inline FloatType GetAssignment(int indiv, int chromosome, int locus, int cluster) const
	{
#ifdef USE_VECTOR
		return assignments[indiv][chromosome][locus][cluster];
#else
		return assignments[GetIdx(indiv, chromosome, locus, cluster)];
#endif
	}

	inline void SetAssignment(int indiv, int chromosome, int locus, int cluster, FloatType val)
	{
#ifdef USE_VECTOR
		assignments[indiv][chromosome][locus][cluster] = val;
#else
		assignments[GetIdx(indiv, chromosome, locus, cluster)] = val;
#endif
	}

	inline int GetNumIndivs() const { return num_indivs; }
	inline int GetNumLoci() const { return num_loci; }
	inline int GetNumClusters() const { return num_clusters; }

	inline int GetIdx(int indiv, int chrom, int locus, int cluster) const
	{
		// Indiv_n   <   l_0:[a: Z_1 .. Z_K | b: Z_1 .. Z_K]      ...      l_L:[a: Z_1 .. Z_K | b: Z_1 .. Z_K]   >
		const int INDIV_START = indiv * GetNumLoci() * GetNumClusters() * NUM_CHROMOSOMES;
		const int LOCUS_START = locus * GetNumClusters() * NUM_CHROMOSOMES;
		const int CLUSTER_START = chrom * GetNumClusters();
		return INDIV_START + LOCUS_START + CLUSTER_START + cluster;
	}

	inline void Update(const BitGenosMatrix& genos, const struct P& p, const struct Q& q);

	inline void Normalize()
	{
		for (int n = 0; n < GetNumIndivs(); ++n) {
			for (int l = 0; l < GetNumLoci(); ++l) {
				FloatType sm_c0 = static_cast<FloatType>(0.0);
				FloatType sm_c1 = static_cast<FloatType>(0.0);
				for (int k = 0; k < GetNumClusters(); ++k) {
					sm_c0 += GetAssignment(n, 0, l, k);
					sm_c1 += GetAssignment(n, 1, l, k);
				}
				for (int k = 0; k < GetNumClusters(); ++k) {
					const FloatType norm_c0 = GetAssignment(n, 0, l, k) / sm_c0;
					const FloatType norm_c1 = GetAssignment(n, 1, l, k) / sm_c1;
					SetAssignment(n, 0, l, k, norm_c0);
					SetAssignment(n, 1, l, k, norm_c1);
				}
			}
		}
	}

	inline double GetEntropy(int indiv, int chromosome, int locus) const
	{
		double ent = -LogFact(GetNumIndivs());

		double sm = 0.0;
		for (int i = 0; i < GetNumClusters(); ++i) {
			const FloatType z_i = GetAssignment(indiv, chromosome, locus, i);
			sm += z_i * log(z_i);
		}
		ent -= GetNumIndivs() * sm;

		sm = 0.0;
		for (int i = 0; i < GetNumClusters(); ++i) {
			for (int x_i = 0; x_i < GetNumClusters(); ++x_i) {
				const FloatType p_i = GetAssignment(indiv, chromosome, locus, i);
				sm += Comb(GetNumClusters(), x_i) * pow(p_i, x_i) * pow(1.0 - p_i, GetNumClusters() - x_i) * LogFact(x_i);
			}
		}
		ent += sm;
		return ent;
	}

	int num_indivs;
	int num_loci;
	int num_clusters;

#ifdef USE_VECTOR
	std::vector<std::vector<std::vector<std::vector<FloatType>>>> assignments;
#else
	FloatType* assignments;
#endif
};

struct P
{
	static constexpr int NUM_PARAMS = 2;

	P() : num_loci(0), num_clusters(0), beta(0), gamma(0)
#ifndef USE_VECTOR
		, freqs(nullptr)
#endif
	{}

	~P()
	{
#ifndef USE_VECTOR
		delete[] freqs; freqs = nullptr;
#endif
	}

	void Init(int num_loci, int num_clusters)
	{
		this->num_loci = num_loci;
		this->num_clusters = num_clusters;
		beta = gamma = static_cast<FloatType>(0.5);

#ifdef USE_VECTOR
		freqs.resize(GetNumClusters());
		for (int k = 0; k < GetNumClusters(); ++k)
			freqs[k].resize(GetNumLoci());
#else
		const int NUM_ELEMS = NUM_PARAMS * GetNumClusters() * GetNumLoci();
		freqs = new FloatType[NUM_ELEMS];
#endif

		// Initialize uniform.
		for (int l = 0; l < GetNumLoci(); ++l)
			for (int k = 0; k < GetNumClusters(); ++k) {
				SetFreq(l, k, 0, static_cast<FloatType>(1.0));
				SetFreq(l, k, 1, static_cast<FloatType>(1.0));
			}
	}

	inline FloatType GetFreq(int num_loci, int num_cluster, int param_idx) const
	{
#ifdef USE_VECTOR
		const auto& p = freqs[num_cluster][num_loci];
		return param_idx == 0 ? p.first : p.second;
#else
		return freqs[GetIdx(num_loci, num_cluster, param_idx)];
#endif
	}

	inline void SetFreq(int num_loci, int num_cluster, int param_idx, FloatType val)
	{
#ifdef USE_VECTOR
		auto& p = freqs[num_cluster][num_loci];
		if (param_idx == 0)
			p.first = val;
		else
			p.second = val;
#else
		freqs[GetIdx(num_loci, num_cluster, param_idx)] = val;
#endif
	}

	inline int GetNumLoci() const { return num_loci; }
	inline int GetNumClusters() const { return num_clusters; }

	inline int GetIdx(int num_loci, int num_cluster, int param_idx) const
	{
		// P   <   l_0:[u: P_1 .. P_K | v: P_1 .. P_K]      ...      l_L[u: P_1 .. P_K | v: P_1 .. P_K]   >
		const int LOCI_START = num_loci * NUM_PARAMS * GetNumClusters();
		const int PARAM_START = param_idx * NUM_PARAMS;
		return LOCI_START + PARAM_START + num_cluster;
	}

	inline void Update(const BitGenosMatrix& genos, const Z& z)
	{
		for (int k = 0; k < GetNumClusters(); ++k) {
			for (int l = 0; l < GetNumLoci(); ++l) {
				FloatType sm_za = static_cast<FloatType>(0.0);
				FloatType sm_zb = static_cast<FloatType>(0.0);
				for (int n = 0; n < genos.GetNumIndivs(); ++n) {
					const int G = genos.GetGeno(n, l);
					const FloatType z_a = z.GetAssignment(n, 0, l, k);
					const FloatType z_b = z.GetAssignment(n, 1, l, k);
					const FloatType z_ab = z_a + z_b;
					sm_za += (G == 1 ? z_a : 0) + (G == 2 ? z_ab : 0);
					sm_zb += (G == 1 ? z_b : 0) + (G == 0 ? z_ab : 0);
				}
				SetFreq(l, k, 0, beta + sm_za);
				SetFreq(l, k, 1, gamma + sm_zb);
			}
		}
	}

	inline void Update2(const BitGenosMatrix& genos, const Z& z)
	{
		for (int l = 0; l < GetNumLoci(); ++l) {
			for (int k = 0; k < GetNumClusters(); ++k) {
				FloatType sm_za = static_cast<FloatType>(0.0);
				FloatType sm_zb = static_cast<FloatType>(0.0);
				for (int n = 0; n < genos.GetNumIndivs(); ++n) {
					const int G = genos.GetGeno(n, l);
					const FloatType z_a = z.GetAssignment(n, 0, l, k);
					const FloatType z_b = z.GetAssignment(n, 1, l, k);
					const FloatType z_ab = z_a + z_b;
					sm_za += (G == 1 ? z_a : 0) + (G == 2 ? z_ab : 0);
					sm_zb += (G == 1 ? z_b : 0) + (G == 0 ? z_ab : 0);
				}
				SetFreq(l, k, 0, beta + sm_za);
				SetFreq(l, k, 1, gamma + sm_zb);
			}
		}
	}

	int num_loci;
	int num_clusters;
	FloatType beta;
	FloatType gamma;

#ifdef USE_VECTOR
	std::vector<std::vector<std::pair<FloatType, FloatType>>> freqs;
#else
	FloatType* freqs;
#endif
};

struct Q
{
	Q() : num_indivs(0), num_clusters(0), alpha(0)
#ifndef USE_VECTOR
		, props(nullptr)
#endif
	{}

	~Q()
	{
#ifndef USE_VECTOR
		delete[] props; props = nullptr;
#endif
		num_indivs = num_clusters = 0;
	}

	inline void Init(int num_indivs, int num_clusters)
	{
		this->num_indivs = num_indivs;
		this->num_clusters = num_clusters;
		alpha = static_cast<FloatType>(1.0 / GetNumClusters());

#ifdef USE_VECTOR
		props.resize(GetNumIndivs());
		for (int n = 0; n < GetNumIndivs(); ++n)
			props[n].resize(GetNumClusters());
#else
		const int NUM_ELEMS = GetNumIndivs() * GetNumClusters();
		props = new FloatType[NUM_ELEMS];
#endif

		// Initialize uniform.
		for (int n = 0; n < GetNumIndivs(); ++n)
			for (int k = 0; k < GetNumClusters(); ++k)
				SetAdmixProp(n, k, static_cast<FloatType>(1.0 / GetNumClusters()));
	}

	inline FloatType GetAdmixProp(int indiv, int cluster) const
	{
#ifdef USE_VECTOR
		return props[indiv][cluster];
#else
		return props[GetIdx(indiv, cluster)];
#endif
	}

	inline void SetAdmixProp(int indiv, int cluster, FloatType val)
	{
#ifdef USE_VECTOR
		props[indiv][cluster] = val;
#else
		props[GetIdx(indiv, cluster)] = val;
#endif
	}

	inline int GetNumIndivs() const { return num_indivs; }
	inline int GetNumClusters() const { return num_clusters; }
	inline int GetIdx(int indiv, int cluster) const { return GetNumClusters() * indiv + cluster; }

	inline void Update(const Z& z)
	{
		for (int n = 0; n < GetNumIndivs(); ++n) {
			for (int k = 0; k < GetNumClusters(); ++k) {
				FloatType sm_z_ab = static_cast<FloatType>(0.0);
				for (int l = 0; l < z.GetNumLoci(); ++l)
					sm_z_ab += z.GetAssignment(n, 0, l, k) + z.GetAssignment(n, 1, l, k);
				SetAdmixProp(n, k, alpha + sm_z_ab);
			}
		}
	}

	inline FloatType GetQ0(int n) const
	{
		FloatType sm = static_cast<FloatType>(0.0);
		for (int k = 0; k < GetNumClusters(); ++k)
			sm += GetAdmixProp(n, k);
		return sm;
	}

	inline int GetIndivCluster(int n) const
	{
		const FloatType Q_0 = GetQ0(n);
		FloatType max_prop = static_cast<FloatType>(-1.0);
		int max_k = -1;
		for (int k = 0; k < GetNumClusters(); ++k) {
			const FloatType PROP = GetAdmixProp(n, k) / Q_0;
			if (max_prop < PROP) {
				max_prop = PROP;
				max_k = k;
			}
		}
		return max_k;
	}

	int num_indivs;
	int num_clusters;
	FloatType alpha;

#ifdef USE_VECTOR
	std::vector<std::vector<FloatType>> props;
#else
	FloatType* props;
#endif
};



void Z::Update(const BitGenosMatrix& genos, const P& p, const Q& q)
{
	for (int n = 0; n < GetNumIndivs(); ++n) {
		for (int l = 0; l < GetNumLoci(); ++l) {
			FloatType q_0 = q.GetQ0(n);
			for (int k = 0; k < GetNumClusters(); ++k) {
				const int G = genos.GetGeno(n, l);
				const FloatType p_u = p.GetFreq(l, k, 0);
				const FloatType p_v = p.GetFreq(l, k, 1);
				const FloatType q_n = q.GetAdmixProp(n, k);

				const FloatType psi_a = (G == 0 ? digamma(p_v) : 0) +
										(G == 1 ? digamma(p_u) : 0) +
										(G == 2 ? digamma(p_u) : 0);

				const FloatType psi_b = (G == 0 ? digamma(p_v) : 0) +
										(G == 1 ? digamma(p_v) : 0) +
										(G == 2 ? digamma(p_u) : 0);

				const FloatType dg_p_uv = digamma(p_u + p_v);
				const FloatType dg_q_n = digamma(q_n);
				const FloatType dg_q_0 = digamma(q_0);
				const FloatType z_a = exp(psi_a - dg_p_uv + dg_q_n - dg_q_0);
				const FloatType z_b = exp(psi_b - dg_p_uv + dg_q_n - dg_q_0);
				SetAssignment(n, 0, l, k, z_a);
				SetAssignment(n, 1, l, k, z_b);
			}
		}
	}
	Normalize();
}



static void ReadGenotypes(std::ifstream& geno_file, uint32_t num_indivs, uint32_t num_loci, BitGenosMatrix& genos)
{
	int bit_num = 0, buff;
	geno_file.read(reinterpret_cast<char*>(&buff), 1);
	logger << Time << " Reading genotypes . . ." << std::endl;
	for (unsigned i = 0; i < num_indivs; ++i) {
		//logger << "#" << i << "     ";
		for (unsigned l = 0; l < num_loci; ++l) {
			if (bit_num == CHAR_BIT) {
				geno_file.read(reinterpret_cast<char*>(&buff), 1);
				bit_num = 0;
			}
			const int G = (buff >> bit_num) & BitGenosMatrix::GENOTYPE_MASK;
			//logger << G << ' ';
			genos.SetGeno(i, l, G);
			bit_num += BitGenosMatrix::BITS_PER_GENOTYPE;
		}
		//logger << std::endl;
	}
	logger << Time << " Reading is done!" << std::endl;
}

static double CalculateLLBO(const BitGenosMatrix& genos, const Z& z, const Q& q, const P& p)
{
	const double LOG_BETA_B_G = LogBeta(p.beta, p.gamma);

	double LLBO = 0.0;
	for (int l = 0; l < genos.GetNumLoci(); ++l) {
		for (int n = 0; n < genos.GetNumIndivs(); ++n) {
			const int G = genos.GetGeno(n, l);
			for (int k = 0; k < z.GetNumClusters(); ++k) {
				// E [[ Z_{nlk}^a ]]      E [[ Z_{nlk}^b ]]
				const FloatType exp_za = genos.GetNumIndivs() * z.GetAssignment(n, 0, l, k);
				const FloatType exp_zb = genos.GetNumIndivs() * z.GetAssignment(n, 1, l, k);
				const double ent_za = z.GetEntropy(n, 0, l);
				const double ent_zb = z.GetEntropy(n, 1, l);

				// E [[ log(1 - P_{lk}) ]]     E [[ log(P_{lk}) ]]
				const FloatType p_u = p.GetFreq(l, k, 0);
				const FloatType p_v = p.GetFreq(l, k, 1);
				const FloatType dg_p_uv = digamma(p_u + p_v);
				const FloatType exp_1_Plk = digamma(p_v) - dg_p_uv;
				const FloatType exp_Plk = digamma(p_u) - dg_p_uv;

				// E [[ log Q_{nk} ]]
				const FloatType exp_Qnk = digamma(q.GetAdmixProp(n, k)) - digamma(q.GetQ0(n));

				const double z_ab = exp_za + exp_zb;
				if (G == 0)
					LLBO += z_ab * (exp_1_Plk);
				else if (G == 2)
					LLBO += z_ab * exp_Plk;
				LLBO += z_ab * exp_Qnk;

				if (G == 1)
					LLBO += exp_za * exp_Plk + exp_zb * exp_1_Plk - ent_za - ent_zb;

				const double BETA_RATIO = LogBeta(p_u, p_v) - LOG_BETA_B_G;
				LLBO += BETA_RATIO + (p.beta - p_u) * exp_Plk + (p.gamma - p_v) * exp_1_Plk;
			}
		}
	}

	const double log_dg_alpha_0 = lgamma(q.alpha * q.GetNumClusters());
	for (int n = 0; n < genos.GetNumIndivs(); ++n) {
		const double q_n0 = q.GetQ0(n);
		LLBO += lgamma(q_n0) - log_dg_alpha_0;
	}
	return LLBO;
}

inline static bool IsConverged(double new_LLBO, double old_LLBO)
{
	//if (new_LLBO < old_LLBO) return true;
	const double DIFF = new_LLBO - old_LLBO;
	return DIFF < LLBO_EPSILON;
}

inline double CalcLogProb(int indiv, int cluster, const BitGenosMatrix& genos, const Z& z, const P& p)
{
	const double LOG_2 = log(2);
	double prob = 0.0;
	for (int l = 0; l < z.GetNumLoci(); ++l) {
		const int G = genos.GetGeno(indiv, l);
		const FloatType p_u = p.GetFreq(l, cluster, 0);
		const FloatType p_lk = p_u / (p_u + p.GetFreq(l, cluster, 1));
		if (G == 0)
			prob += log(1.0 - p_lk) + log(1.0 - p_lk);
		else if (G == 1)
			prob += LOG_2 + log(p_lk) + log(1.0 - p_lk);
		else
			prob += log(p_lk) + log(p_lk);
	}
	return prob;
}

inline static void DumpVarParams(const BitGenosMatrix& genos, const Z& z, const Q& q, const P& p)
{
	logger << Time << " Dumping variational parameters . . ." << std::endl;

	std::ofstream prop_file(DUMP_PATH + "props.txt");
	if (!prop_file.is_open()) {
		logger << Time << ' ' << warning << " Could not open output file for dumping variational parameters!" << std::endl;
		return;
	}

	for (int n = 0; n < q.GetNumIndivs(); ++n) {
		prop_file << n << "     ";
		const FloatType Q_0 = q.GetQ0(n);
		for (int k = 0; k < q.GetNumClusters(); ++k)
			prop_file << q.GetAdmixProp(n, k) / Q_0 << ' ';
		prop_file << "      Q: " << q.GetIndivCluster(n) << "     LogProbs: ";

		double max_f = -std::numeric_limits<double>::max(); int max_k = -100;
		for (int k = 0; k < z.GetNumClusters(); ++k) {
			const double LOG_PROB = CalcLogProb(n, k, genos, z, p);
			prop_file << LOG_PROB << ' ';
			if (max_f < LOG_PROB) {
				max_f = LOG_PROB;
				max_k = k;
			}
		}
		prop_file << "     Z:" << max_k << std::endl;
	}

	logger << Time << " Dumping is done!" << std::endl;
}

inline static bool DumpFreqs(std::string path, FreqsVector& freqs)
{
	std::ofstream freqs_file(path);
	if (!freqs_file.is_open())
		return false;

	// Dump configs.
	freqs_file << "NUM_CLUSTERS: " << freqs.size() << std::endl;
	freqs_file << "NUM_LOCI: " << freqs[0].size() << std::endl;
	freqs_file << std::endl;

	// Dump frequencies.
	for (unsigned k = 0; k < freqs.size(); ++k) {
		for (unsigned l = 0; l < freqs[0].size(); ++l)
			freqs_file << freqs[k][l].first << ' ';
		freqs_file << std::endl;
		for (unsigned l = 0; l < freqs[0].size(); ++l)
			freqs_file << (freqs[k][l].second - freqs[k][l].first) << ' ';
		freqs_file << std::endl;
		for (unsigned l = 0; l < freqs[0].size(); ++l)
			freqs_file << (1.0 - freqs[k][l].second) << ' ';
		freqs_file << std::endl << std::endl;
	}
	freqs_file.close();
	return true;
}

inline static bool ReadFreqsFromTextFile(std::string path, FreqsVector& freqs)
{
	freqs.clear();

	std::ifstream freqs_file(path);
	if (!freqs_file.is_open())
		return false;

	// Read configs.
	int num_clusters, num_loci;
	std::string tmp_str;
	freqs_file >> tmp_str >> num_clusters >> tmp_str >> num_loci;

	// Read frequencies.
	freqs.resize(num_clusters);
	for (int k = 0; k < num_clusters; ++k) {
		freqs[k].resize(num_loci);
		double freq;
		for (int l = 0; l < num_loci; ++l) {
			freqs_file >> freq;
			freqs[k][l].first = freq;
		}
		for (int l = 0; l < num_loci; ++l) {
			freqs_file >> freq;
			freqs[k][l].second = freqs[k][l].first + freq;
		}
		for (int l = 0; l < num_loci; ++l)
			freqs_file >> freq;			// Skip unused frequencies.
	}
	return true;
}

static bool IsFileExist(const std::string& path)
{
#if defined __linux__ || defined __CYGWIN__
	return std::filesystem::exists(path);
#else
	return std::experimental::filesystem::exists(path);
#endif
}

inline static bool InitFreqs(FreqsVector& freqs, int num_clusters, int num_loci)
{
	if (IsFileExist(FREQS_PATH)) {
		logger << Time << " Freqs file exists . . ." << std::endl;
		return ReadFreqsFromTextFile(FREQS_PATH, freqs);
	}

	logger << Time << " Making random frequencies . . ." << std::endl;

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	freqs.resize(num_clusters);
	for (int k = 0; k < static_cast<int>(freqs.size()); ++k) {
		freqs[k].resize(num_loci);
		for (int l = 0; l < num_loci; ++l) {
			double p0 = unif(rd);
			double p1 = unif(rd);
			if (p0 > p1)
				std::swap(p0, p1);
			freqs[k][l].first = p0;
			freqs[k][l].second = p1;
		}
	}
	return DumpFreqs(FREQS_PATH, freqs);
}

inline bool InitGenos(BitGenosMatrix& genos, const FreqsVector& freqs, int num_indivs, int num_loci)
{
	if (IsFileExist(GENOS_PATH)) {
		logger << Time << " Genotype file exists . . ." << std::endl;
		return genos.ReadFromTextFile(GENOS_PATH);
	}

	logger << Time << " Making genotype file . . ." << std::endl;

	const int NUM_CLUSTERS = static_cast<int>(freqs.size());
	genos.Init(num_indivs * NUM_CLUSTERS, num_loci, NUM_CLUSTERS);

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	for (int k = 0; k < NUM_CLUSTERS; ++k)
		for (int i = 0; i < num_indivs; ++i)
			for (int l = 0; l < num_loci; ++l) {
				const double R = unif(rd);
				const int G = (R < freqs[k][l].first) ? 0 : (R < freqs[k][l].second) ? 1 : 2;
				const int INDIV_IDX = num_indivs * k + i;
				genos.SetGeno(INDIV_IDX, l, G);
			}
	return genos.DumpText(GENOS_PATH);
}

inline static void LogIterations(int itr, int MAX_ITERS, double old_llbo)
{
	(void) MAX_ITERS;
#ifdef USE_LLBO
	//logger << "LLBO:" << old_LLBO << std::endl;
#endif

	if (itr % ITER_REPORT == 0)
		logger << Time << " itr #" << itr << "     LLBO:" <<
#ifdef USE_LLBO
		old_llbo
#else
		"NOT USED"
#endif
		<< std::endl;
}

inline static void CalcAcc(const Q& q)
{
	std::vector<std::vector<int>> cluster_count(q.GetNumClusters(),
		std::vector<int>(q.GetNumClusters(), 0));

	// Count clusters.
	const int NUM_INDIVS_IN_CLUSTER = q.GetNumIndivs() / q.GetNumClusters();
	for (int k = 0; k < q.GetNumClusters(); ++k)
		for (int i = 0; i < NUM_INDIVS_IN_CLUSTER; ++i)
			++cluster_count[k][q.GetIndivCluster(k * NUM_INDIVS_IN_CLUSTER + i)];

	std::vector<int> perm(q.GetNumClusters());
	for (int i = 0; i < q.GetNumClusters(); ++i)
		perm[i] = i;

	// Calculate accuracies and find maximum one for report.
	std::vector<double> accs;
	do {
		double sm = 0.0;
		for (int k = 0; k < q.GetNumClusters(); ++k)
			sm += cluster_count[k][perm[k]];
		const double ACC = sm / static_cast<double>(q.GetNumIndivs());
		accs.push_back(ACC);
	} while (std::next_permutation(perm.begin(), perm.end()));

	logger << Time << " Acc: " << *std::max_element(accs.begin(), accs.end()) << "   [ ";
	for (const auto& a : accs)
		logger << a << ' ';
	logger << ']' << std::endl;
}



int main(int argc, char** argv)
{
	logger << std::endl << std::endl;
	logger << "===== VB STRUCT =====" << std::endl;
	logger << "Start : " << Time << std::endl;
	logger << "MAX_ITERS    : " << MAX_ITERS << std::endl;
	logger << "MIN_ITERS    : " << MIN_ITERS << std::endl;
	//logger << "LLBO_UPDATE  : " << LLBO_UPDATE << std::endl;
	logger << "LLBO_EPSILON : " << LLBO_EPSILON << std::endl;
	logger << "sizeof(int)  : " << sizeof(int) << std::endl;
	logger << "USE VECTOR   : " << (USE_VECTOR ? "TRUE" : "FALSE") << std::endl;
	logger << std::endl;

	if (argc < 2) {
		logger << "File name is required!" << std::endl;
		return 1;
	}

#ifdef READ_GENOTYPES_FROM_BINARY_FILE
	std::ifstream geno_file(argv[1], std::ios_base::binary);
	if (!geno_file.is_open()) {
		logger << "Could not open input file!" << std::endl;
		return 2;
	}

	// Read header.
	logger << "Genotype file `" << argv[1] << "' is opened!" << std::endl;
	uint32_t num_clusters, num_loci, num_indivs;
	geno_file.read(reinterpret_cast<char*>(&num_clusters), sizeof(num_clusters));
	geno_file.read(reinterpret_cast<char*>(&num_loci), sizeof(num_loci));
	geno_file.read(reinterpret_cast<char*>(&num_indivs), sizeof(num_indivs));
	logger << "Header ->     <K:" << num_clusters << ">    <L:" << num_loci << ">   <N:" << num_indivs << ">" << std::endl;

	// Read genotypes.
	BitGenosMatrix genos(num_indivs, num_loci);
	ReadGenotypes(geno_file, num_indivs, num_loci, genos);
#elifdef MAKE_RANDOM_FREQS
	//Make frequencies.
	FreqsVector freqs;
	InitFreqs(freqs, NUM_CLUSTERS, NUM_LOCI);

	// Make genotypes.
	InitGenos(genos, freqs, NUM_INDIVS, NUM_LOCI);
#else
	logger << Time << " Reading from text file -> " << argv[1] << std::endl;
	BitGenosMatrix genos;
	if (!genos.ReadFromTextFile(argv[1])) {
		logger << "Could not open `" << argv[1] << "' file!" << std::endl;
		return 2;
	}
#endif

	logger << Time << "Reading genotype file is done!" << std::endl;
	logger << "  NumIndivs:   " << genos.GetNumIndivs() << std::endl;
	logger << "  NumLoci:     " << genos.GetNumLoci() << std::endl;
	logger << "  NumClusters: " << genos.GetNumClusters() << std::endl;

	// Initialize parameters.
	logger << Time << " Initialize P, Z, and Q . . ." << std::endl;
	P p; p.Init(genos.GetNumLoci(), genos.GetNumClusters());
	Z z; z.Init(genos.GetNumIndivs(), genos.GetNumLoci(), genos.GetNumClusters());
	Q q; q.Init(genos.GetNumIndivs(), genos.GetNumClusters());

#ifdef MAKE_RANDOM_FREQS
	freqs.clear();		// Clear useless frequencies.
#endif
	
	logger << Time << " Start iterations . . ." << std::endl;

	double old_LLBO =
#ifdef USE_LLBO
		CalculateLLBO(genos, z, q, p);
#else
		0;
#endif

	for (int itr = 0; itr < MAX_ITERS; ++itr) {
		LogIterations(itr, MAX_ITERS, old_LLBO);

		p.Update(genos, z);			// Update P
		z.Update(genos, p, q);		// Update Z
		q.Update(z);				// Update Q

#ifdef USE_LLBO
		const double NEW_LLBO = CalculateLLBO(genos, z, q, p);
		if (IsConverged(NEW_LLBO, old_LLBO) && itr > MIN_ITERS) {
			logger << Time << " Converged at #" << itr << " iteration!         "
				<< "new LLBO:" << NEW_LLBO << "     old LLBO:" << old_LLBO << std::endl;
			break;
		}
		old_LLBO = NEW_LLBO;
#endif
	}

	DumpVarParams(genos, z, q, p);		// Dump variational parameters.
	CalcAcc(q);							// Report accuracy of clustering.
	logger << "End : " << Time << std::endl << std::endl;
	return 0;
}

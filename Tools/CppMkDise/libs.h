#ifndef LIBS_H_
#define LIBS_H_

#include <algorithm>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>



#define WRITE_ALL_TO_FILE	1



/// Individulas and population
typedef std::vector<std::vector<double>> SubFreqs;
typedef std::vector<SubFreqs> FreqsVect;
typedef std::vector<std::vector<int>> IndivGenotype;

// Disease model
struct SNPRec
{
	int locus;		/// locus in { 0, ..., NUM_LOCI }
	int SNP;		/// SNP in   { 0, 1, 2 }
};

typedef std::vector<SNPRec> SNPRecVect;
typedef std::vector<SNPRecVect> DisModel;
typedef std::vector<std::vector<double>> Params;



constexpr int NUM_CHROMOSOMES = 2;
constexpr int REPORT_INDIV = 500;

static int NUM_INDIVS;
static int NUM_LOCI;
static int NUM_CLUSTERS;
static int NUM_DIFFS;
static std::string FMT;
static std::string FREQS_PATH;
static std::string GENOS_PATH;
static std::string LABEL_PATH;



static double GetMultiVarSigmoid(const std::vector<double>& params, const std::vector<int>& xs)
{
	double sm = params[0];		// First parameter is bias.
	for (int i = 0; i < static_cast<int>(xs.size()); ++i)
		sm += params[i + 1] * xs[i];
	return 1.0 / (1.0 + exp(-sm));
}

static std::string GetCTime()
{
	time_t tm;
	time(&tm);
	std::string tm_str(ctime(&tm));
	char c = tm_str[tm_str.size() - 1];
	while (c == '\n' || c == '\r') {
		tm_str = tm_str.substr(0, tm_str.size() - 1);
		c = tm_str[tm_str.size() - 1];
	}
	return tm_str + ' ';
}

static void MakeFreqs(FreqsVect& freqs)
{
	std::cout << GetCTime() << "MakeFreqs() -> " << FREQS_PATH << std::endl;

	std::vector<int> diffs;
	std::random_device rd;
	std::uniform_real_distribution<double> unif_diff(0, NUM_LOCI);
	for (int d = 0; d < NUM_DIFFS; ++d) {
		int l = int(unif_diff(rd));
		while (std::find(diffs.begin(), diffs.end(), l) != diffs.end())
			l = int(unif_diff(rd));
		diffs.push_back(l);
	}
	std::sort(diffs.begin(), diffs.end());

	std::uniform_real_distribution<double> unif(0, 1);
	freqs.resize(NUM_CLUSTERS);
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		freqs[k].resize(NUM_CHROMOSOMES);
		for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
			unsigned d_idx = 0;
			freqs[k][c].resize(NUM_LOCI);
			for (int l = 0; l < NUM_LOCI; ++l) {
				if (k == 0)
					freqs[k][c][l] = unif(rd);
				else {
					if (d_idx < diffs.size() && l == diffs[d_idx]) {
						freqs[k][c][l] = unif(rd);
						++d_idx;
					} else
						freqs[k][c][l] = freqs[0][c][l];
				}
			}
		}
	}
	std::cout << "freqs.size()    " << freqs.size() << std::endl;
	std::cout << "freqs[0].size() " << freqs[0].size() << std::endl;
	std::cout << "freqs[0][0].size() " << freqs[0][0].size() << std::endl;

#ifdef WRITE_ALL_TO_FILE
	std::ofstream freqs_file(FREQS_PATH);
	if (!freqs_file.is_open()) {
		std::cout << "Could not open freqs file! -> `" << FREQS_PATH << '\'' << std::endl;
		return;
	}

	freqs_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	freqs_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	freqs_file << std::endl;

	for (const auto& d : diffs)
		freqs_file << d << ' ';
	freqs_file << std::endl << std::endl;

	for (const auto& sub : freqs) {
		for (const auto& ch : sub) {
			for (double p : ch)
				freqs_file << p << ' ';
			freqs_file << std::endl;
		}
		freqs_file << std::endl;
	}
#endif
}

static bool IsInDiseaseModels(int l, const std::vector<DisModel>& models)
{
	for (const DisModel& dis : models)
		for (const SNPRecVect& vect : dis)
			for (const SNPRec& rec : vect)
				if (rec.locus == l)
					return true;

	return false;
}

static void MakeDiseaseFreqs(FreqsVect& freqs, const std::vector<DisModel>& models)
{
	std::cout << GetCTime() << "MakeFreqs() -> " << FREQS_PATH << std::endl;

	std::vector<int> diffs;
	std::random_device rd;
	std::uniform_real_distribution<double> unif_diff(0, NUM_LOCI);
	for (int d = 0; d < NUM_DIFFS; ++d) {
		int l = int(unif_diff(rd));
		while (std::find(diffs.begin(), diffs.end(), l) != diffs.end() || IsInDiseaseModels(l, models))
			l = int(unif_diff(rd));
		diffs.push_back(l);
	}
	std::sort(diffs.begin(), diffs.end());

	std::uniform_real_distribution<double> unif(0, 1);
	freqs.resize(NUM_CLUSTERS);
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		freqs[k].resize(NUM_CHROMOSOMES);
		for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
			unsigned d_idx = 0;
			freqs[k][c].resize(NUM_LOCI);
			for (int l = 0; l < NUM_LOCI; ++l) {
				if (k == 0)
					freqs[k][c][l] = unif(rd);
				else {
					if (d_idx < diffs.size() && l == diffs[d_idx]) {
						freqs[k][c][l] = unif(rd);
						++d_idx;
					} else
						freqs[k][c][l] = freqs[0][c][l];
				}
			}
		}
	}
	std::cout << "freqs.size()    " << freqs.size() << std::endl;
	std::cout << "freqs[0].size() " << freqs[0].size() << std::endl;
	std::cout << "freqs[0][0].size() " << freqs[0][0].size() << std::endl;

#ifdef WRITE_ALL_TO_FILE
	std::ofstream freqs_file(FREQS_PATH);
	if (!freqs_file.is_open()) {
		std::cout << "Could not open freqs file! -> `" << FREQS_PATH << '\'' << std::endl;
		return;
	}

	freqs_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	freqs_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	freqs_file << std::endl;

	for (const auto& d : diffs)
		freqs_file << d << ' ';
	freqs_file << std::endl << std::endl;

	for (const auto& sub : freqs) {
		for (const auto& ch : sub) {
			for (double p : ch)
				freqs_file << p << ' ';
			freqs_file << std::endl;
		}
		freqs_file << std::endl;
	}
#endif
}

static void MakeGenos(const FreqsVect& freqs)
{
	std::cout << GetCTime() << "MakeGenos()" << std::endl;

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	std::ofstream genos_file(GENOS_PATH);

#ifdef WRITE_ALL_TO_FILE
	std::ofstream my_genos_file(GENOS_PATH + "-faststr.txt");
	my_genos_file << "NUM_INDIVS: " << NUM_INDIVS << std::endl;
	my_genos_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	my_genos_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	my_genos_file << std::endl;

	std::vector<int> genos(NUM_LOCI);
#endif

	const int TOTAL_INDIVS = NUM_CLUSTERS * NUM_INDIVS;
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		for (int i = 0; i < NUM_INDIVS; ++i) {
			const int IND_IDX = k * NUM_INDIVS + i;
			if (IND_IDX % REPORT_INDIV == 0) {
				const int PERCENT = static_cast<int>((IND_IDX + 1) / static_cast<float>(TOTAL_INDIVS) * 100.0);
				std::cout << GetCTime() << "IND_IDX:" << IND_IDX + 1 << "    " << PERCENT << '%' << std::endl;
			}
			for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
				genos_file << "IND_" << IND_IDX + 1 << ' ';
				for (int l = 0; l < NUM_LOCI; ++l) {
					const double P = freqs[k][c][l];
					const double R = unif(rd);
					const int G = R < P ? 0 : 1;
					genos_file << G;
					if (l + 1 < NUM_LOCI)
						genos_file << ' ';

#ifdef WRITE_ALL_TO_FILE
					if (c == 0)
						genos[l] = G;
					else
						genos[l] += G;
#endif
				}
				if (IND_IDX + 1 < TOTAL_INDIVS || c + 1 < NUM_CHROMOSOMES)
					genos_file << '\n';
			}

#ifdef WRITE_ALL_TO_FILE
			for (int l = 0; l < NUM_LOCI; ++l)
				my_genos_file << genos[l];
			if (IND_IDX + 1 < TOTAL_INDIVS)
				my_genos_file << std::endl;
#endif
		}
	}
}

static bool IsCaseIndiv(const IndivGenotype& dis_genos, const DisModel& model, const std::vector<double>& params)
{
	std::vector<int> xs;
	for (int i = 0; i < static_cast<int>(model.size()); ++i) {
		xs.push_back(1);
		const SNPRecVect& vect = model[i];
		for (int r = 0; r < static_cast<int>(vect.size()); ++r) {
			const int L = vect[r].locus;
			const int G = dis_genos[0][L] + dis_genos[1][L];
			xs[i] *= (G == vect[r].SNP) ? 1 : 0;
		}
	}
	const double RES = GetMultiVarSigmoid(params, xs);
	return RES > 0.5;
}



static void MakeDiseGenos1(const FreqsVect& freqs, const std::vector<DisModel>& models, const Params& params)
{
	std::cout << GetCTime() << "MakeGenos()" << std::endl;

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	std::ofstream genos_file(GENOS_PATH);
	std::ofstream label_file(LABEL_PATH);

	IndivGenotype dise_genos(2);
	for (int c = 0; c < NUM_CHROMOSOMES; ++c)
		dise_genos[c].resize(NUM_LOCI);

#ifdef WRITE_ALL_TO_FILE
	std::ofstream my_genos_file(GENOS_PATH + "-faststr.txt");
	my_genos_file << "NUM_INDIVS: " << NUM_INDIVS << std::endl;
	my_genos_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	my_genos_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	my_genos_file << std::endl;

	std::vector<int> genos(NUM_LOCI);
#endif

	const int TOTAL_INDIVS = NUM_CLUSTERS * NUM_INDIVS;
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		for (int i = 0; i < NUM_INDIVS; ++i) {
			bool is_case;
			const int IND_IDX = k * NUM_INDIVS + i;
			if (IND_IDX % REPORT_INDIV == 0) {
				const int PERCENT = static_cast<int>((IND_IDX + 1) / static_cast<float>(TOTAL_INDIVS) * 100.0);
				std::cout << GetCTime() << "IND_IDX:" << IND_IDX + 1 << "    " << PERCENT << '%' << std::endl;
			}
			while (true) {
				for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
					for (int l = 0; l < NUM_LOCI; ++l) {
						const double P = freqs[k][c][l];
						const double R = unif(rd);
						const int G = R < P ? 0 : 1;
						dise_genos[c][l] = G;

#ifdef WRITE_ALL_TO_FILE
						if (c == 0)
							genos[l] = G;
						else
							genos[l] += G;
#endif
					}
				}

				is_case = IsCaseIndiv(dise_genos, models[k], params[k]);
				if ((i < NUM_INDIVS / 2)) {
					if (!is_case)
						continue;
				} else if (is_case)
					continue;
				break;
			}

			if (is_case)
				label_file << "YES" << std::endl;
			else
				label_file << "NO" << std::endl;

			for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
				genos_file << "IND_" << IND_IDX + 1 << ' ';
				for (int l = 0; l < NUM_LOCI; ++l) {
					genos_file << dise_genos[c][l];
					if (l + 1 < NUM_LOCI)
						genos_file << ' ';
				}
				if (IND_IDX + 1 < TOTAL_INDIVS || c + 1 < NUM_CHROMOSOMES)
					genos_file << '\n';
			}
#ifdef WRITE_ALL_TO_FILE
			for (int l = 0; l < NUM_LOCI; ++l)
				my_genos_file << genos[l];
			if (IND_IDX + 1 < TOTAL_INDIVS)
				my_genos_file << std::endl;
#endif
		}
	}
}



static void MakeIndivGenos(const SubFreqs& freqs, IndivGenotype& indiv)
{
	if (static_cast<int>(indiv.size()) != NUM_CHROMOSOMES || static_cast<int>(indiv[0].size()) != NUM_LOCI) {
		std::cout << "Indiv size is not compatible!" << std::endl;
		exit(1);
	}

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
		for (int l = 0; l < static_cast<int>(freqs[0].size()); ++l) {
			const double P = freqs[c][l];
			const double R = unif(rd);
			const int G = R < P ? 0 : 1;
			indiv[c][l] = G;
		}
	}
}

static void WriteStruct(int ind_idx, int total_indivs, const IndivGenotype& indiv, std::ofstream& genos_file)
{
	for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
		genos_file << "IND_" << ind_idx + 1 << ' ';
		for (int l = 0; l < NUM_LOCI; ++l) {
			genos_file << indiv[c][l];
			if (l + 1 < NUM_LOCI)
				genos_file << ' ';
		}
		if (ind_idx + 1 < total_indivs || c + 1 < NUM_CHROMOSOMES)
			genos_file << '\n';
	}
}

static void WriteMyStruct(int ind_idx, int total_indivs, const IndivGenotype& indiv, std::ofstream& my_genos_file)
{
#ifdef WRITE_ALL_TO_FILE
	for (int l = 0; l < NUM_LOCI; ++l)
		my_genos_file << (indiv[0][l] + indiv[1][l]);
	if (ind_idx + 1 < total_indivs)
		my_genos_file << std::endl;
#endif
}

static void MakeDiseGenos2(const FreqsVect& freqs, const std::vector<DisModel>& models, const Params& params)
{
	std::cout << GetCTime() << "MakeGenos()" << std::endl;

	std::ofstream genos_file(GENOS_PATH);
	std::ofstream label_file(LABEL_PATH);

	IndivGenotype dise_genos(2);
	for (int c = 0; c < NUM_CHROMOSOMES; ++c)
		dise_genos[c].resize(NUM_LOCI);

#ifdef WRITE_ALL_TO_FILE
	std::ofstream my_genos_file(GENOS_PATH + "-faststr.txt");
	my_genos_file << "NUM_INDIVS: " << NUM_INDIVS << std::endl;
	my_genos_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	my_genos_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	my_genos_file << std::endl;
#endif

	const int TOTAL_INDIVS = NUM_CLUSTERS * NUM_INDIVS;
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		std::vector<IndivGenotype> cntr_indivs;
		const int NUM_DISEASE_INDIVS = NUM_INDIVS / 2;
		int i = 0;
		while (i < NUM_DISEASE_INDIVS) {
			MakeIndivGenos(freqs[k], dise_genos);
			if (!IsCaseIndiv(dise_genos, models[k], params[k])) {
				if (static_cast<int>(cntr_indivs.size()) < NUM_INDIVS - NUM_DISEASE_INDIVS)
					cntr_indivs.push_back(dise_genos);
				continue;
			}

			const int INDIV_IDX = k * NUM_INDIVS + i;
			WriteStruct(INDIV_IDX, TOTAL_INDIVS, dise_genos, genos_file);
			WriteMyStruct(INDIV_IDX, TOTAL_INDIVS, dise_genos, my_genos_file);
			label_file << "YES" << std::endl;
			++i;
		}

		for (int j = 0; j < static_cast<int>(cntr_indivs.size()); ++j) {
			const int INDIV_IDX = k * NUM_INDIVS + i;
			WriteStruct(INDIV_IDX, TOTAL_INDIVS, cntr_indivs[j], genos_file);
			WriteMyStruct(INDIV_IDX, TOTAL_INDIVS, cntr_indivs[j], my_genos_file);
			label_file << "NO" << std::endl;
			++i;
		}

		while (i < NUM_INDIVS) {
			MakeIndivGenos(freqs[k], dise_genos);
			if (IsCaseIndiv(dise_genos, models[k], params[k]))
				continue;

			const int INDIV_IDX = k * NUM_INDIVS + i;
			WriteStruct(i, TOTAL_INDIVS, dise_genos, genos_file);
			WriteMyStruct(i, TOTAL_INDIVS, dise_genos, my_genos_file);
			label_file << "NO" << std::endl;
			++i;
		}
	}
}



static void AddDisease(IndivGenotype& dise_genos,
#if defined WRITE_ALL_TO_FILE
		std::vector<int>& genos,
#endif
		const DisModel& model)
{
	for (const SNPRecVect& vect : model)
		for (const SNPRec& rec : vect) {
			if (rec.SNP == 0) {
				dise_genos[0][rec.locus] = dise_genos[1][rec.locus] = 0;
			} else if (rec.SNP == 2) {
				dise_genos[0][rec.locus] = dise_genos[1][rec.locus] = 1;
			} else {
				dise_genos[0][rec.locus] = 0;
				dise_genos[1][rec.locus] = 1;
			}

#if defined WRITE_ALL_TO_FILE
			genos[rec.locus] = rec.SNP;
#endif
		}
}

static void MakeDiseGenos(const FreqsVect& freqs, const std::vector<DisModel>& models, const Params& params)
{
	std::cout << GetCTime() << "MakeGenos()" << std::endl;

	std::random_device rd;
	std::uniform_real_distribution<double> unif(0, 1);
	std::ofstream genos_file(GENOS_PATH);
	std::ofstream label_file(LABEL_PATH);

	IndivGenotype dise_genos(2);
	for (int c = 0; c < NUM_CHROMOSOMES; ++c)
		dise_genos[c].resize(NUM_LOCI);

#ifdef WRITE_ALL_TO_FILE
	std::ofstream my_genos_file(GENOS_PATH + "-faststr.txt");
	my_genos_file << "NUM_INDIVS: " << NUM_INDIVS << std::endl;
	my_genos_file << "NUM_LOCI: " << NUM_LOCI << std::endl;
	my_genos_file << "NUM_CLUSTERS: " << NUM_CLUSTERS << std::endl;
	my_genos_file << std::endl;

	std::vector<int> genos(NUM_LOCI);
#endif

	const int TOTAL_INDIVS = NUM_CLUSTERS * NUM_INDIVS;
	for (int k = 0; k < NUM_CLUSTERS; ++k) {
		for (int i = 0; i < NUM_INDIVS; ++i) {
			bool is_case;
			const int IND_IDX = k * NUM_INDIVS + i;
			if (IND_IDX % REPORT_INDIV == 0) {
				const int PERCENT = static_cast<int>((IND_IDX + 1) / static_cast<float>(TOTAL_INDIVS) * 100.0);
				std::cout << GetCTime() << "IND_IDX:" << IND_IDX + 1 << "    " << PERCENT << '%' << std::endl;
			}
			while (true) {
				for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
					for (int l = 0; l < NUM_LOCI; ++l) {
						const double P = freqs[k][c][l];
						const double R = unif(rd);
						const int G = R < P ? 0 : 1;
						dise_genos[c][l] = G;

#ifdef WRITE_ALL_TO_FILE
						if (c == 0)
							genos[l] = G;
						else
							genos[l] += G;
#endif
					}
				}

				if (i < NUM_INDIVS / 2)
					AddDisease(dise_genos, genos, models[k]);

				is_case = IsCaseIndiv(dise_genos, models[k], params[k]);
				if (i < NUM_INDIVS / 2) {
					if (!is_case)
						continue;
				} else if (is_case)
					continue;
				break;
			}

			if (is_case)
				label_file << "YES" << std::endl;
			else
				label_file << "NO" << std::endl;

			for (int c = 0; c < NUM_CHROMOSOMES; ++c) {
				genos_file << "IND_" << IND_IDX + 1 << ' ';
				for (int l = 0; l < NUM_LOCI; ++l) {
					genos_file << dise_genos[c][l];
					if (l + 1 < NUM_LOCI)
						genos_file << ' ';
				}
				if (IND_IDX + 1 < TOTAL_INDIVS || c + 1 < NUM_CHROMOSOMES)
					genos_file << '\n';
			}
#ifdef WRITE_ALL_TO_FILE
			for (int l = 0; l < NUM_LOCI; ++l)
				my_genos_file << genos[l];
			if (IND_IDX + 1 < TOTAL_INDIVS)
				my_genos_file << std::endl;
#endif
		}
	}
}

#endif


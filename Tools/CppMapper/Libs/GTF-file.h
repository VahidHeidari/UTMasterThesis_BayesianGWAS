#ifndef GTF_FILE_H_
#define GTF_FILE_H_

#include <cctype>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <vector>

class GTFFile
{
public:
	typedef std::vector<std::string> GTFRecs;

	enum GTF_FILED {
		GTF_FIELD_SEQ_ID,
		GTF_FIELD_SOURCE,
		GTF_FIELD_TYPE,
		GTF_FIELD_START,
		GTF_FIELD_END,
		GTF_FIELD_SCORE,
		GTF_FIELD_STRAND,
		GTF_FIELD_PHASE,
		GTF_FIELD_ATTRIBUTES,
	};

	GTFFile(std::string input_file) : input_file_path(input_file) {}

	bool Read(std::string out_file_path, const std::set<std::string>& gene_names)
	{
		std::ifstream in_file(input_file_path);
		if (!in_file.is_open()) {
			std::cout << "Could not open input file! -> " << input_file_path << std::endl;
			return false;
		}

		std::ofstream out_file(out_file_path);
		if (!out_file.is_open()) {
			std::cout << "Could not open output file! -> " << out_file_path << std::endl;
			return false;
		}

		std::string line;
		std::vector<std::string> sp;
		while (std::getline(in_file, line)) {
			if (!line.size() || line[0] == '#')
				continue;

			SplitLine('\t', line, sp);
			if (sp[GTF_FIELD_TYPE] != "gene")
				continue;

			const std::string& GENE_NAME = GetGeneID(sp);
			if (gene_names.find(GENE_NAME) != gene_names.end())
				genes.push_back(sp);
		}
		std::cout << "Num Genes:" << genes.size() << std::endl;
		return true;
	}

//private:
	void SplitLine(char separator, const std::string& line, GTFRecs& sp) const
	{
		sp.clear();
		size_t st = 0;
		while (st < line.size()) {
			size_t ed = line.find(separator, st);
			if (ed == std::string::npos)
				ed = line.size();
			sp.push_back(Trim(line.substr(st, ed - st)));
			st = ed + 1;
		}
	}

	std::string Trim(const std::string& s) const
	{
		size_t st = 0;
		while (st < s.size() && isspace(s[st]))
			++st;
		if (st == s.size())
			return "";

		size_t ed = s.size() - 1;
		while (ed > 0 && isspace(s[ed]))
			--ed;

		return s.substr(st, ed - st + 1);
	}

	std::string GetGeneID(const GTFRecs& sp) const
	{
		GTFRecs attr_sp, attrs;
		SplitLine(';', sp[GTF_FIELD_ATTRIBUTES], attr_sp);
		for (unsigned i = 0; i < sp.size(); ++i) {
			SplitLine(' ', sp[i], attrs);
			if (attrs[0] == "gene_id")
				return attrs[1].substr(1, attrs[1].size() - 3);
		}

		return "";
	}

	std::string input_file_path;
	std::vector<GTFRecs> genes;
};

#endif


#ifndef SAM_FILE_H_
#define SAM_FILE_H_

#include <cstdlib>
#include <cctype>
#include <vector>



typedef std::vector<std::string> SAMRec;



enum SAM_FIELDS {
	QNAME = 0,
	FLAG = 1,
	RNAME = 2,
	POS = 3,
	MAPQ = 4,
	CIGAR = 5,
	RNEXT = 6,
	PNEXT = 7,
	TLEN = 8,
	SEQ = 9,
	QUAL = 10
};

enum SAM_FLAGS {
	FLAGS_MULTI_SEGMENTS = 0x0001,
	FLAGS_ALIGNED = 0x0002,
	FLAGS_UNMAPPED = 0x0004,
	FLAGS_NEXT_SEGMENT_UNMAPPED = 0x0008,
	FLAGS_REVERSED_COMPLEMENTED = 0x0010,
	FLAGS_NEXT_SEGMENT_REVERSED = 0x0020,
	FLAGS_FIRST_SEGMENT = 0x0040,
	FLAGS_LAST_SEGMENT = 0x0080,
	FLAGS_SECONDARY_ALIGNMENT = 0x0100,
	FLAGS_NOT_PASSING_QUALITY_CONTROL = 0x0200,
	FLAGS_PCR_OR_OPTICAL_DUPLICATED = 0x0400
};



struct CIGAROps
{
	std::vector<char> ops;
	std::vector<int> ops_len;
};



static SAMRec Split(const std::string& line)
{
	SAMRec res;
	size_t tab_pos = 0;
	while (true) {
		const size_t cur_tab_pos = line.find('\t', tab_pos);
		if (cur_tab_pos == std::string::npos)
			break;

		res.push_back(line.substr(tab_pos, cur_tab_pos - tab_pos));
		tab_pos = cur_tab_pos + 1;
	}

	res.push_back(line.substr(tab_pos, line.size()));
	return res;
}

static inline int GetFlag(const SAMRec& rec)
{
	return atoi(rec[SAM_FIELDS::FLAG].c_str());
}

static inline bool HasFlag(const SAMRec& rec, SAM_FLAGS flag)
{
	return (GetFlag(rec) & flag) == flag;
}

static inline const std::string& GetSeq(const SAMRec& rec)
{
	return rec[SAM_FIELDS::SEQ];
}

static inline const std::string& GetRName(const SAMRec& rec)
{
	return rec[SAM_FIELDS::RNAME];
}

static inline int GetPos(const SAMRec& rec)
{
	return atoi(rec[SAM_FIELDS::POS].c_str());
}

static inline const std::string& GetCIGAR(const SAMRec& rec)
{
	return rec[SAM_FIELDS::CIGAR];
}

static inline void SplitCIGARFromStr(const std::string& cigar, CIGAROps& ops)
{
	ops.ops.clear();
	ops.ops_len.clear();
	size_t i = 0;
	while (i < cigar.size()) {
		size_t j = i;
		while (j < cigar.size() && isdigit(cigar[j]))
			++j;

		ops.ops_len.push_back(atoi(cigar.substr(i, j).c_str()));
		ops.ops.push_back(cigar[j]);
		i = j + 1;
	}
}

static inline void SplitCIGAR(const SAMRec& rec, CIGAROps& ops)
{
	const std::string& cigar = GetCIGAR(rec);
	SplitCIGARFromStr(cigar, ops);
}

#endif

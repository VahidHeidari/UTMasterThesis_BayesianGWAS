#ifndef BAM_FILE_H_
#define BAM_FILE_H_

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gzip-stream.h"

static constexpr double MAX_PROB_WRONG_POS = 0.01;

struct RefInfo
{
	RefInfo() {}
	RefInfo(const RefInfo& o) : l_name(o.l_name), l_ref(o.l_ref), name(o.name) {}

	unsigned int GetSize() const { return 8 + name.size(); }

	std::ostream& Print(std::ostream& o=std::cout) const
	{
		o << "l_name   : " << l_name << std::endl;
		o << "name     : " << name   << std::endl;
		o << "l_ref    : " << l_ref  << std::endl;
		return o;
	}

	uint32_t l_name;
	uint32_t l_ref;
	std::string name;
};

struct AuxiliaryRec
{
};

struct AlignmentRec
{
	static constexpr int HDR_SIZE = 36;

	/// CIGAR operations.
	enum CIGAR_OP {
		CIGAR_OP_MATCH = 0,			/// M
		CIGAR_OP_INSERT,			/// I
		CIGAR_OP_SKIP_RIGION,		/// N
		CIGAR_OP_SOFT_CLIP,			/// S
		CIGAR_OP_HARD_CLIP,			/// H
		CIGAR_OP_PADDING,			/// P
		CIGAR_OP_SEQ_MATCH,			/// =
		CIGAR_OP_SEQ_MISMATCH,		/// X
	};

	AlignmentRec() {}
	AlignmentRec(const char* blk) { Read(blk); }
	AlignmentRec(const AlignmentRec& o)
		: block_size(o.block_size), ref_ID(o.ref_ID), pos(o.pos),
		l_read_name(o.l_read_name), mapq(o.mapq), bin(o.bin),
		n_cigar_op(o.n_cigar_op), flag(o.flag), l_seq(o.l_seq),
		next_ref_ID(o.next_ref_ID), next_pos(o.next_pos), tlen(o.tlen),
		read_name(o.read_name), cigar(o.cigar), seq(o.seq), qual(o.qual)
		//auxs(o.auxs)
	{}

	void Read(const char* buff)
	{
		CastT(block_size,   0, buff);
		CastT(ref_ID,       4, buff);
		CastT(pos,          8, buff);
		CastT(l_read_name, 12, buff);
		CastT(mapq,        13, buff);
		CastT(bin,         14, buff);
		CastT(n_cigar_op,  16, buff);
		CastT(flag,        18, buff);
		CastT(l_seq,       20, buff);
		CastT(next_ref_ID, 24, buff);
		CastT(next_pos,    28, buff);
		CastT(tlen,        32, buff);
		//ReadString(f, read_name, l_read_name);
		//ReadString(f, cigar, 4 * n_cigar_op);
		//ReadString(f, seq, (l_seq + 1) / 2);
		//ReadString(f, qual, l_seq);

		// DOTO: Skip Auxiliary records for now!
		//int rem_bytes = block_size + 4 - 36 - l_read_name - (4 * n_cigar_op) - ((l_seq + 1) / 2) - l_seq;
		//f.seekg(rem_bytes, std::ios_base::cur);
	}

	uint32_t GetCIGARBin(int idx) const { return *reinterpret_cast<const uint32_t*>(&cigar[idx * 4]); }
	int GetCIGARLength(int idx) const { return static_cast<int>(GetCIGARBin(idx) >> 4 & 0x0FFFFFFF); }
	int GetCIGAROp(int idx) const { return static_cast<int>(GetCIGARBin(idx) & 0x0F); }
	unsigned int GetSize() const { return block_size + 4; }
	bool IsUnmapped() const { return (flag & 0x4) || (flag & 0x8); }
	bool IsSecondaryAlignment() const { return flag & 0x800; }

	bool HasSoftOrHardClipping() const
	{
		for (unsigned i = 0; i < n_cigar_op; ++i) {
			const int OP = GetCIGAROp(i);
			if (OP == CIGAR_OP::CIGAR_OP_HARD_CLIP || OP == CIGAR_OP::CIGAR_OP_SOFT_CLIP)
				return true;
		}

		return false;
	}

	std::string GetCIGARStr() const
	{
		std::string cigar_str;
		for (unsigned i = 0; i < n_cigar_op; ++i)
			cigar_str += std::to_string(GetCIGARLength(i)) + "MIDNSHP=X"[GetCIGAROp(i)];
		return cigar_str;
	}

	std::string GetSeqStr() const
	{
		std::string seq_str;
		for (unsigned i = 0; i < l_seq; ++i)
			seq_str += "=ACMGRSVTWYHKDBN"[(i % 2 ? seq[i / 2] : seq[i / 2] >> 4) & 0x0F];
		return seq_str;
	}

	bool IsOmittedQual() const
	{
		for (unsigned i = 0; i < l_seq; ++i) {
			if (qual[i] != static_cast<char>(0xFF))
				return false;
		}

		return true;
	}

	std::string GetQualStr() const
	{
		if (IsOmittedQual())
			return "*";

		std::string qual_str;
		qual_str.resize(l_seq);
		for (unsigned i = 0; i < l_seq; ++i)
			qual_str[i] = qual[i] + 33;
		return qual_str;
	}

	std::string GetNextRefID() const
	{
		if (next_ref_ID == -1)
			return "*";

		if (next_ref_ID == ref_ID)
			return "=";

		return std::to_string(next_ref_ID + 1);
	}

	int GetPosWithSoftClip() const
	{
		if (n_cigar_op == 0)
			return pos;

		if (GetCIGAROp(0) == CIGAR_OP::CIGAR_OP_SOFT_CLIP)
			return pos - GetCIGARLength(0);

		return pos;
	}

	int GetSeqLength() const
	{
		if (n_cigar_op == 0)
			return l_seq;

		int l = l_seq;
		for (int i = 0; i < n_cigar_op; ++i) {
			const int OP = GetCIGAROp(i);
			if (OP == CIGAR_OP::CIGAR_OP_SKIP_RIGION)
				l += GetCIGARLength(i);
			else if (OP == CIGAR_OP::CIGAR_OP_INSERT)
				l -= GetCIGARLength(i);
		}
		return l;
	}

	void Print(std::ostream& o=std::cout, std::string idnt="") const
	{
		o << idnt << "block_size    : " << block_size    << std::endl;
		o << idnt << "ref_ID        : " << ref_ID        << std::endl;
		o << idnt << "pos           : " << pos           << std::endl;
		o << idnt << "l_read_name   : " << static_cast<int>(l_read_name) << std::endl;
		o << idnt << "mapq          : " << static_cast<int>(mapq) << std::endl;
		o << idnt << "bin           : " << bin           << std::endl;
		o << idnt << "n_cigar_op    : " << n_cigar_op    << std::endl;
		o << idnt << "flag          : " << flag          << std::endl;
		o << idnt << "l_seq         : " << l_seq         << std::endl;
		o << idnt << "next_ref_ID   : " << next_ref_ID   << std::endl;
		o << idnt << "next_pos      : " << next_pos      << std::endl;
		o << idnt << "tlen          : " << tlen          << std::endl;
		o << idnt << "read_name     : " << read_name     << std::endl;
		o << idnt << "cigar         : " << GetCIGARStr() << std::endl;
		o << idnt << "seq           : " << GetSeqStr()   << std::endl;
		o << idnt << "qual          : " << GetQualStr()  << std::endl;
	}

	uint32_t block_size;
	int32_t ref_ID;
	int32_t pos;
	uint8_t l_read_name;
	uint8_t mapq;
	uint16_t bin;
	uint16_t n_cigar_op;
	uint16_t flag;
	uint32_t l_seq;
	int32_t next_ref_ID;
	int32_t next_pos;
	int32_t tlen;
	std::string read_name;
	std::string cigar;
	std::string seq;
	std::string qual;
	//std::vector<AuxiliaryRec> auxs;		/// TODO: Implement auxiliary tags!
};

struct Region
{
	Region() : start_pos(0), end_pos(0) {}
	Region(int st, int ed) : start_pos(st), end_pos(ed) {}
	Region(const Region& o) : start_pos(o.start_pos), end_pos(o.end_pos) {}

	bool HasIntersection(int start_pos, int end_pos) const
	{
		return start_pos <= this->end_pos && this->start_pos <= end_pos;
	}

	int start_pos;
	int end_pos;
};

typedef std::vector<Region> RegionsVect;
typedef std::map<std::string, RegionsVect> RegionsMap;

class BAMFile
{
public:
	BAMFile(std::string file_path) : l_text(0), n_ref(0), alg_cur_idx(0), gzip_stream(file_path) { memset(magic, 0, 4); }

	void Read()
	{
		ReadHeader();				// Read header.
		ReadReferenceInfos();		// Read reference informations.
		ReadAlignments();			// Read alignments.
	}

	void Read(const RegionsMap& regions)
	{
		ReadHeader();				// Read header.
		ReadReferenceInfos();		// Read reference informations.
		ReadAlignments(regions);	// Read alignments.
	}

	bool DumpSAM(std::string out_path) const
	{
		std::ofstream o(out_path);
		if (!o.is_open()) {
			std::cout << "Could not open `" << out_path << "' for dump SAM!" << std::endl;
			return false;
		}

		// Dump header.
		o << "@HD\tVN:1.3\tSO:coordinate" << std::endl;
		for (const auto& r : ref_infos)
			o << "@SQ\tSN:" << r.name.substr(0, r.name.size() - 1) << "\tLN:" << r.l_ref << std::endl;
		o << "@CO\tCppGeneExpr generated output" << std::endl;

		// Dump alignment.
		for (const auto& a : alignments) {
			const RefInfo& r = ref_infos[a.ref_ID];
			o << a.read_name.substr(0, a.read_name.size() - 1) << '\t'
				<< a.flag << '\t'
				<< r.name.substr(0, r.name.size() - 1) << '\t'
				<< a.pos + 1 << '\t'
				<< static_cast<int>(a.mapq) << '\t'
				<< a.GetCIGARStr() << '\t'
				<< a.GetNextRefID() << '\t'
				<< a.next_pos + 1 << '\t'
				<< a.tlen << '\t'
				<< a.GetSeqStr() << '\t'
				<< a.GetQualStr() << '\t'
				<< std::endl;
		}
		return true;
	}

	std::ostream& Print(std::ostream& o=std::cout) const
	{
		std::string magic_str(magic, 4);
		o << "Header:" << std::endl;
		o << "----------" << std::endl;
		o << "magic   : " << magic_str << std::endl;
		o << "l_text  : " << l_text    << std::endl;
		o << "text    : " << text      << std::endl;
		o << "n_ref   : " << n_ref     << std::endl;

		o << std::endl;
		o << "Ref infos (" << ref_infos.size() << ") :" << std::endl;
		o << "----------" << std::endl;
		for (const auto& r : ref_infos) {
			r.Print(o);
			o << "----------" << std::endl;
		}

		o << std::endl;
		o << "Alignments (" << alignments.size() << ") :" << std::endl;
		o << "--------------------" << std::endl;
		for (const auto& a : alignments) {
			a.Print(o, "  ");
			o << "--------------------" << std::endl;
		}
		return o;
	}

	void ReadHeaderAndReferenceInfos()
	{
		ReadHeader();
		ReadReferenceInfos();
	}

	bool GetNextAlignmentRec(AlignmentRec& alg)
	{
		if (!alg_cur_idx)
			alg_cur_idx = GetStartOfAlignmentIdx();

		const char* blk = gzip_stream.GetBlock(alg_cur_idx, AlignmentRec::HDR_SIZE);
		if (gzip_stream.IsEOF())
			return false;

		alg.Read(blk);
		ReadString(alg.read_name, alg_cur_idx + AlignmentRec::HDR_SIZE, alg.l_read_name);
		ReadString(alg.cigar,     alg_cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name, 4 * alg.n_cigar_op);
		ReadString(alg.seq,       alg_cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op), (alg.l_seq + 1) / 2);
		ReadString(alg.qual,      alg_cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op) + ((alg.l_seq + 1) / 2), alg.l_seq);
		alg_cur_idx += alg.GetSize();
		return true;
	}

	const std::string& GetRefName(const AlignmentRec& alg) const { return ref_infos[alg.ref_ID].name; }

	std::string GetRefNameStr(const AlignmentRec& alg) const
	{
		const std::string& ref_name_z = GetRefName(alg);
		const std::string ref_name = ref_name_z.substr(0, ref_name_z.size() - 1);
		return ref_name;
	}

private:
	void ReadHeader()
	{
		memcpy(magic, gzip_stream.GetBlock(0, 4), 4);
		Read(l_text, 4);
		ReadString(text, 8, l_text);
		Read(n_ref, 8 + l_text);
	}

	void ReadReferenceInfos()
	{
		int cur_idx = GetHeaderSize();
		for (unsigned i = 0; i < n_ref; ++i) {
			RefInfo ref;
			Read(ref.l_name, cur_idx);
			ReadString(ref.name, cur_idx + 4, ref.l_name);
			Read(ref.l_ref, cur_idx + 4 + ref.l_name);
			ref_infos.push_back(ref);
			cur_idx += ref.GetSize();
		}
	}

	void ReadAlignments()
	{
		int num_alignments = 0;
		int cur_idx = GetStartOfAlignmentIdx();
		while (true) {
			const char* blk = gzip_stream.GetBlock(cur_idx, AlignmentRec::HDR_SIZE);
			if (gzip_stream.IsEOF())
				break;

			AlignmentRec alg(blk);
			ReadString(alg.read_name, cur_idx + AlignmentRec::HDR_SIZE, alg.l_read_name);
			ReadString(alg.cigar,     cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name, 4 * alg.n_cigar_op);
			ReadString(alg.seq,       cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op), (alg.l_seq + 1) / 2);
			ReadString(alg.qual,      cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op) + ((alg.l_seq + 1) / 2), alg.l_seq);
			alignments.push_back(alg);
			cur_idx += alg.GetSize();
			++num_alignments;
		}
		std::cout << "NumAlignments: " << num_alignments << std::endl;
	}

	void ReadAlignments(const RegionsMap& regions)
	{
		unsigned int num_alignments = 0;
		unsigned int cur_idx = GetStartOfAlignmentIdx();
		while (true) {
			const char* blk = gzip_stream.GetBlock(cur_idx, AlignmentRec::HDR_SIZE);
			if (gzip_stream.IsEOF())
				break;

			// Deserializing alignment.
			AlignmentRec alg(blk);
			ReadString(alg.read_name, cur_idx + AlignmentRec::HDR_SIZE, alg.l_read_name);
			ReadString(alg.cigar,     cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name, 4 * alg.n_cigar_op);
			ReadString(alg.seq,       cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op), (alg.l_seq + 1) / 2);
			ReadString(alg.qual,      cur_idx + AlignmentRec::HDR_SIZE + alg.l_read_name + (4 * alg.n_cigar_op) + ((alg.l_seq + 1) / 2), alg.l_seq);

			// Check filtering conditions.
			if (!IsFiltered(alg, regions)) {
				alignments.push_back(alg);
				++num_alignments;
			}
			cur_idx += alg.GetSize();		// Advance current index.
		}
		std::cout << "NumAlignments: " << num_alignments << std::endl;
	}

	inline bool IsFiltered(const AlignmentRec& alg, const RegionsMap& regions) const
	{
		// Filter unmapped or secondary alignments.
		if (alg.IsUnmapped() || alg.IsSecondaryAlignment())
			return true;

		// Filter soft clip or hard clip alignments.
		if (alg.HasSoftOrHardClipping())
			return true;

		// Filter low quality mappings.
		const double PROB_WRONG_POS = pow(10.0, -alg.mapq / 10.0);
		if (PROB_WRONG_POS > MAX_PROB_WRONG_POS)
			return true;

		// Filter out of bound alignments.
		if (!HasIntersection(alg, regions))
			return true;

		return false;
	}

	template <class T>
	void Read(T& t, int start_idx)
	{
		t = *reinterpret_cast<const T*>(gzip_stream.GetBlock(start_idx, sizeof(T)));
	}

	void ReadString(std::string& out_str, int start_idx, int str_length)
	{
		out_str = std::string(gzip_stream.GetBlock(start_idx, str_length), str_length);
	}

	unsigned int GetHeaderSize() const { return 12 + l_text; }
	unsigned int GetStartOfAlignmentIdx() const { return GetHeaderSize() + GetRefInfosSize(); }

	unsigned int GetRefInfosSize() const
	{
		unsigned int ref_infos_size = 0;
		for (const auto& r : ref_infos)
			ref_infos_size += r.GetSize();
		return ref_infos_size;
	}

	bool HasIntersection(const AlignmentRec& alg, const RegionsMap& regions) const
	{
		const auto& v_itr = regions.find(GetRefNameStr(alg));
		if (v_itr == regions.end())
			return false;

		const int ALG_START = alg.GetPosWithSoftClip() + 1;
		const int ALG_END = ALG_START + alg.GetSeqLength();
		const auto& r_itr = v_itr->second;
		for (const auto& r : r_itr)
			if (r.HasIntersection(ALG_START, ALG_END))
				return true;

		return false;
	}

	char magic[4];
	uint32_t l_text;
	std::string text;
	uint32_t n_ref;
	int alg_cur_idx;
	GzipStream gzip_stream;
	std::vector<RefInfo> ref_infos;
	std::vector<AlignmentRec> alignments;
};

#endif


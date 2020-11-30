#ifndef GZIP_BLOCK_H_
#define GZIP_BLOCK_H_

#include <exception>
#include <sstream>
#include <string>
#include <vector>

#include <zlib.h>

#include "util.h"

class GzipBlock
{
public:
	static constexpr int HDR_SIZE = 12 + 6;
	static constexpr int FOOTER_SIZE = 8;
	static constexpr int EOF_BLOCK_SIZE = HDR_SIZE + 2 + FOOTER_SIZE;

	inline GzipBlock() {}
	inline GzipBlock(std::ifstream& f) { Read(f); }

	void Read(std::ifstream& f)
	{
		// Read header fileds.
		std::string buff;
		buff.resize(HDR_SIZE);
		f.read(&buff[0], HDR_SIZE);
		CastT(ID1,    0, buff.data());
		CastT(ID2,    1, buff.data());
		CastT(CM,     2, buff.data());
		CastT(FLG,    3, buff.data());
		CastT(MTIME,  4, buff.data());
		CastT(XFL,    8, buff.data());
		CastT(OS,     9, buff.data());
		CastT(XLEN,  10, buff.data());
		CastT(SI1,   12, buff.data());
		CastT(SI2,   13, buff.data());
		CastT(SLEN,  14, buff.data());
		CastT(BSIZE, 16, buff.data());

		// Read compressed data.
		const int COMP_SIZE = BSIZE - XLEN - 19;
		const int HDR_COMP_FOOTER_SIZE = HDR_SIZE + COMP_SIZE + FOOTER_SIZE;
		CDATA.resize(HDR_COMP_FOOTER_SIZE);
		memcpy(&CDATA[0], buff.data(), HDR_SIZE);
		f.read(reinterpret_cast<char*>(&CDATA[HDR_SIZE]), COMP_SIZE + FOOTER_SIZE);

		// Read footer fields.
		const int FOOTER_OFF = HDR_SIZE + COMP_SIZE;
		CastT(CRC32, FOOTER_OFF, CDATA.data());
		CastT(ISIZE, FOOTER_OFF + 4, CDATA.data());
		//std::cout << "s:" << CDATA.size() << "    c:" << CDATA.capacity() << std::endl;
	}

	inline bool IsEOF() const
	{
		return ID1 == 0x1f && ID2 == 0x8b && CM == 8 && FLG == 4 &&
			MTIME == 0 && XFL == 0 && OS == 0xff && XLEN == 6 &&
			SI1 == 66 && SI2 == 67 && SLEN == 2 && BSIZE == 0x001b &&
			CRC32 == 0 && ISIZE == 0 &&
			CDATA.size() == EOF_BLOCK_SIZE && CDATA[HDR_SIZE] == '\x03' && CDATA[HDR_SIZE + 1] == '\0';
	}

	inline std::ostream& Print(std::ostream& os = std::cout) const
	{
		os << "ID1    : " << std::hex << static_cast<int>(ID1) << std::endl;
		os << "ID2    : " << std::hex << static_cast<int>(ID2) << std::endl;
		os << "CM     : " << std::hex << static_cast<int>(CM) << std::endl;
		os << "FLG    : " << std::hex << static_cast<int>(FLG) << std::endl;
		os << "MTIME  : " << std::dec << MTIME << std::endl;
		os << "XFL    : " << std::hex << static_cast<int>(XFL) << std::endl;
		os << "OS     : " << std::hex << static_cast<int>(OS) << std::endl;
		os << "XLEN   : " << std::dec << XLEN << std::endl;
		os << "SI1    : " << std::hex << static_cast<int>(SI1) << std::endl;
		os << "SI2    : " << std::hex << static_cast<int>(SI2) << std::endl;
		os << "SLEN   : " << std::dec << SLEN << std::endl;
		os << "BSIZE  : " << std::dec << BSIZE << std::endl;
		os << "CDATA  : " << std::dec << CDATA.size() << std::endl;
		os << "CRC32  : " << std::hex << CRC32 << std::endl;
		os << "ISIZE  : " << std::dec << ISIZE << std::endl;
		return os;
	}

	bool Uncompress(std::string& uncomp_buff)
	{
		z_stream zs;                        // z_stream is zlib's control structure
		memset(&zs, 0, sizeof(zs));

		static constexpr int MOD_GZIP_ZLIB_WINDOWSIZE = 15;
		if (inflateInit2(&zs, MOD_GZIP_ZLIB_WINDOWSIZE + 16) != Z_OK) {
			std::cout << "inflateInit failed while decompressing." << std::endl;
			return false;
		}

		zs.next_in = (Bytef*)CDATA.data();
		zs.avail_in = CDATA.size();

		// get the decompressed bytes blockwise using repeated calls to inflate
		uncomp_buff.resize(ISIZE);
		zs.next_out = reinterpret_cast<Bytef*>(&uncomp_buff[0]);
		zs.avail_out = ISIZE;
		int ret = inflate(&zs, 0);
		inflateEnd(&zs);

		if (ret != Z_STREAM_END) {          // an error occurred that was not EOF
			std::cout << "Exception during zlib decompression: (" << ret << ") " << zs.msg << std::endl;
			return false;
		}

		return true;
	}

	uint8_t ID1;
	uint8_t ID2;
	uint8_t CM;
	uint8_t FLG;
	uint32_t MTIME;
	uint8_t XFL;
	uint8_t OS;
	uint16_t XLEN;
	uint8_t SI1;
	uint8_t SI2;
	uint16_t SLEN;
	uint16_t BSIZE;
	std::string CDATA;
	uint32_t CRC32;
	uint32_t ISIZE;
};

class GzipStream
{
public:
	GzipStream(std::string file_path) : is_EOF(false), cur_idx(0), file_path(file_path) {}

	const char* GetBlock(unsigned start_idx, unsigned sz)
	{
		// Read first block.
		if (!blk.size()) {
			f.open(file_path, std::ios_base::binary | std::ios_base::in);
			if (!f.is_open()) {
				std::cout << "Could not open input gzip block file!" << std::endl;
				return nullptr;
			}

			gzip_blk.Read(f);
			gzip_blk.Uncompress(blk);
			cur_idx = 0;
			return GetBlock(start_idx, sz);
		}

		// Read new block.
		if (blk.size() + cur_idx == start_idx) {
			gzip_blk.Read(f);
			if (gzip_blk.IsEOF()) {
				is_EOF = true;
				return nullptr;
			}

			gzip_blk.Uncompress(blk);
			cur_idx = start_idx;
			return blk.data();
		}

		// Read next block and append to buffer.
		if (blk.size() + cur_idx < start_idx + sz) {
			std::string buff;
			while (blk.size() + cur_idx < start_idx + sz) {
				gzip_blk.Read(f);
				if (gzip_blk.IsEOF()) {
					is_EOF = true;
					break;
				}
				gzip_blk.Uncompress(buff);
				blk.append(buff.data(), buff.size());
			}
		}

		return blk.data() + start_idx - cur_idx;
	}

	bool IsEOF() const
	{
		return is_EOF;
	}

	unsigned GetCurIdx() const { return cur_idx; }

private:
	bool is_EOF;
	unsigned cur_idx;
	std::string file_path;
	std::ifstream f;
	GzipBlock gzip_blk;
	std::string blk;
};

#endif


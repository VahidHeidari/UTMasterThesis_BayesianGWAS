#include "PackedArray.h"
#include "unit-test.h"

#define WORD_BITS		8		/// Word lenght in bits
#define LENGTH			3		/// Number of words

static inline std::ostream& operator<<(std::ostream& os, const PackedArray& pa)
{
	const int LEN = (int)pa.lengthByte;
	for (int i = 0; i < LEN; ++i) {
		char a = pa.charArray[i];
		os << std::setfill('0') << std::setw(2) << std::hex << ((int)a & 0xFF);
		if (i + 1 < LEN)
			os << ' ';
	}
	return os;
}

void TestPackedArray()
{
	// Test allocation.
	PackedArray pa;
	pa.defineBits(WORD_BITS, LENGTH);
	pa.allocateArray();
	memset(pa.charArray, 0, (int)pa.lengthByte);
	std::cout << pa << std::endl;
	CheckEQ((int)pa.wordLength, WORD_BITS);
	CheckEQ((int)pa.length, LENGTH);
	CheckEQ((int)pa.lengthByte, (int)sizeof(uint) + (LENGTH - 1) * WORD_BITS / 8);

	// Test assignment.
	pa.writePacked(0, 3);			CheckEQ((int)pa[0], 3);
	pa.writePacked(0, 7);			CheckEQ((int)pa[0], 7);
	pa.writePacked(0, 8);			CheckEQ((int)pa[0], 8);
	pa.writePacked(0, 0xff);		CheckEQ((int)pa[0], 0xff);
	std::cout << pa << std::endl;

	pa.writePacked(1, 3);			CheckEQ((int)pa[1], 3);
	pa.writePacked(1, 7);			CheckEQ((int)pa[1], 7);
	pa.writePacked(1, 8);			CheckEQ((int)pa[1], 8);
	pa.writePacked(1, 0xff);		CheckEQ((int)pa[1], 0xff);
	std::cout << pa << std::endl;

	pa.writePacked(2, 3);			CheckEQ((int)pa[2], 3);
	pa.writePacked(2, 7);			CheckEQ((int)pa[2], 7);
	pa.writePacked(2, 8);			CheckEQ((int)pa[2], 8);
	pa.writePacked(2, 0xff);		CheckEQ((int)pa[2], 0xff);
	std::cout << pa << std::endl;

	pa.writePacked(3, 1);
	std::cout << pa << std::endl;

	// Deallocate packed array memory.
	pa.deallocateArray();
}


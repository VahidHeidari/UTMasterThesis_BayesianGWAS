#ifndef UNIT_TEST_H_
#define UNIT_TEST_H_

#include <iostream>

#define TEST(FUNC)	do {											\
	extern void Test##FUNC();										\
	std::cout << "---- START " << #FUNC << "----" << std::endl;		\
	Test##FUNC();													\
	std::cout << "---- END " << #FUNC << "----" << std::endl;		\
} while (false)

template <class L, class R>
bool CheckEQ(const L& l, const R& r)
{
	if (l == r) {
		std::cout << "OK" << std::endl;
		return true;
	}

	std::cout << "Failed!    " << l << " != " << r << std::endl;
	return false;
}

template <class L>
bool CheckTrue(const L& l)
{
	return CheckEQ<L, bool>(l, true);
}

template <class L>
bool CheckFalse(const L& l)
{
	return CheckEQ(l, false);
}

#endif


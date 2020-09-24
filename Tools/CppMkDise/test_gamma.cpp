#include <iostream>

#include "asa147.hpp"

double gammainc(double a, double x)
{
	int e;
	const double RES = gammds(x, a, &e);
	if (e)
		std::cout << "    Error!" << std::endl;
	return RES;
}


double gammaincc(double a, double x)
{
	return 1 - gammainc(a, x);
}

int main()
{
	std::cout << "Test GAMMA" << std::endl;
	int e;
	std::cout << "Gamma(1, 2) = " << gammds(1, 2, &e) << std::endl;
	std::cout << "Gamma(2, 1) = " << gammds(2, 1, &e) << std::endl;
	std::cout << std::endl;

	std::cout << "GammaInc(1, 2) = " << gammainc(1, 2) << std::endl;
	std::cout << "GammaIncc(1, 2) = " << gammaincc(1, 2) << std::endl;
	std::cout << std::endl;

	std::cout << "GammaInc(1, 3) = " << gammainc(1, 3) << std::endl;
	std::cout << "GammaIncc(1, 3) = " << gammaincc(1, 3) << std::endl;
	return 0;
}


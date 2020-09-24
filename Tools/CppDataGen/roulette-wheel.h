#ifndef ROULETTE_WHEEL_H_
#define ROULETTE_WHEEL_H_

#include <random>
#include <vector>

#define SIZE_OF_ARRAY(ARR)		(sizeof(ARR) / sizeof(ARR[0]))



class RouletteWheel
{
public:
	RouletteWheel(std::vector<double> new_chances, bool is_correct_last = false);
	RouletteWheel(const double* new_chances, int chances_size, bool is_correct_last = false);
	RouletteWheel(const RouletteWheel& other);

	RouletteWheel& operator=(const RouletteWheel& other);

	int GetNext();
	void CheckAndInitChances(std::vector<double> new_chances, bool is_correct_last);
	int GetNumChances() const;

	std::uniform_real<double> unif;
	std::vector<double> chances;
};



void TestRouletteWheel();

#endif

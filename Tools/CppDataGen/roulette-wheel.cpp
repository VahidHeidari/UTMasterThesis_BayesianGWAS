#include "roulette-wheel.h"

#include <algorithm>

#include "logger.h"
#include "sim-pars.h"

RouletteWheel::RouletteWheel(std::vector<double> new_chances, bool is_correct_last)
	: unif(0.0, 1.0)
{
	CheckAndInitChances(new_chances, is_correct_last);
}

RouletteWheel::RouletteWheel(const double* new_chances, int chances_size, bool is_correct_last)
	: unif(0.0, 1.0)
{
	std::vector<double> vect_chances;
	for (int i = 0; i < chances_size; ++i)
		vect_chances.push_back(new_chances[i]);

	CheckAndInitChances(vect_chances, is_correct_last);
}

RouletteWheel::RouletteWheel(const RouletteWheel& other)
{
	*this = other;
}

RouletteWheel& RouletteWheel::operator=(const RouletteWheel& other)
{
	unif = other.unif;
	chances = other.chances;
	return *this;
}

int RouletteWheel::GetNext()
{
	std::random_device dev;
	double r = unif(dev);
	double sum = 0;
	for (unsigned i = 0; i < chances.size(); ++i) {
		if (sum <= r && r < sum + chances[i])
			return i;

		sum += chances[i];
	}
	return chances.size() - 1;
}

void RouletteWheel::CheckAndInitChances(std::vector<double> new_chances, bool is_correct_last)
{
	double sum = 0.0;
	for (unsigned i = 0; i < new_chances.size(); ++i)
		sum += new_chances[i];

	double diff = std::abs(1.0 - sum);
	if (diff > EPSILON)
		logger << warning << "Sum of chances in roulette wheel is not equal to one! (sum:"
				<< sum << "    diff:" << diff << ")." << std::endl;

	if (sum > 1.0)
		logger << warning << "Sum of chances is greater than 1!" << std::endl;

	chances = new_chances;
	if (is_correct_last && diff > EPSILON && sum <= 1.0)
		chances.push_back(diff);

	if (sum > 1.0)
		unif = std::uniform_real<double>(0, sum);
}

int RouletteWheel::GetNumChances() const
{
	return chances.size();
}



static void DrawSamples(const double* chances, int chances_size, int num_samples = 10000)
{
	RouletteWheel rw(chances, chances_size, true);
	chances_size = rw.GetNumChances();
	std::vector<int> samples(chances_size);
	for (int i = 0; i < num_samples; ++i) {
		int s = rw.GetNext();
		++samples[s];
	}

	for (int i = 0; i < chances_size; ++i) {
		double ch = (double)samples[i] / (double)(num_samples);
		logger << "Chances " << rw.chances[i] << "     Drawned : " << ch << std::endl;
	}
	logger << std::endl;
}

void TestRouletteWheel()
{
	static const double test_1[] = { 0.3, 0.6, 0.1 };
	static const double test_2[] = { 0.1, 0.3, 0.2, 0.4 };
	static const double test_3[] = { 0.2, 0.3, 0.2, 0.1, 0.2 };
	static const double test_4[] = { 0.2, 0.3, 0, 0.4, 0.1, 0 };
	static const double test_5[] = { 0.7 };
	static const double test_6[] = { 0.3, 0.8 };

	DrawSamples(test_1, SIZE_OF_ARRAY(test_1));
	DrawSamples(test_2, SIZE_OF_ARRAY(test_2));
	DrawSamples(test_3, SIZE_OF_ARRAY(test_3));
	DrawSamples(test_4, SIZE_OF_ARRAY(test_4));
	DrawSamples(test_5, SIZE_OF_ARRAY(test_5));
	DrawSamples(test_6, SIZE_OF_ARRAY(test_6));
}

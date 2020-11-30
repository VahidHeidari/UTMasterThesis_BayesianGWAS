#ifndef LOGGER_H_
#define LOGGER_H_

#include <ctime>
#include <fstream>
#include <iostream>
#include <string>

extern std::ofstream log_file;

struct Log {} extern logger;
struct TimeTag {} extern Time;
struct Warning {} extern warning;

template <typename T>
static Log& operator<<(Log& l, const T& t)
{
	if (log_file.is_open()) {
		log_file << t;
		log_file.flush();
	}

	std::cout << t;
	std::cout.flush();
	return l;
}

static Log& operator<<(Log& l, std::ostream& (*fp)(std::ostream&))
{
	if (log_file.is_open()) {
		fp(log_file);
		log_file.flush();
	}

	fp(std::cout);
	std::cout.flush();
	return l;
}

template <>
Log& operator<< <TimeTag>(Log& l, const TimeTag& /*t*/)
{
	time_t start_time = time(nullptr);

#if defined __linux__ || defined __CYGWIN__
	char* buff = ctime(&start_time);
#else
	enum { TIME_BUFF_SIZE = 30 };
	char buff[TIME_BUFF_SIZE];
	ctime_s(buff, TIME_BUFF_SIZE, &start_time);
#endif

	std::string time_str(buff);
	if (time_str[time_str.size() - 1] == '\n')
		time_str = time_str.substr(0, time_str.size() - 1);
	l << time_str << ' ';
	return l;
}

template <>
Log& operator<< <Warning>(Log& l, const Warning& /*w*/)
{
	l << "WARNING: ";
	return l;
}

#endif

#include "logger.h"

static std::string GetLogFileName()
{
	time_t raw_time;
	time(&raw_time);

	struct tm* time_info = localtime(&raw_time);

	enum { TIME_BUFF_SIZE = 30 };
	char buff[TIME_BUFF_SIZE];
	strftime(buff, TIME_BUFF_SIZE, "log-%Y%m%d.txt", time_info);
	return std::string(buff);
}

std::ofstream log_file(GetLogFileName(), std::ios_base::app | std::ios_base::out);
Log logger;
TimeTag Time;
Warning warning;

#include "logger.h"

#ifdef __linux__
std::ofstream log_file("linux_log_file.txt", std::ios_base::app | std::ios_base::out);
#else
std::ofstream log_file("bin\\log_file.txt", std::ios_base::app | std::ios_base::out);
#endif

Log logger;
TimeTag Time;
Warning warning;

#ifndef UTIL_H_
#define UTIL_H_

#include <iostream>
#include <string>

#include <zlib.h>

template<class T>
static void CastT(T& t, int offset, const char* buff)
{
	t = *reinterpret_cast<const T*>(&buff[offset]);
}

#endif


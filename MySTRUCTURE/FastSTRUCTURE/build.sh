#!/bin/bash

# Build types
#BUILD_TYPE=Debug
BUILD_TYPE=Release
#BUILD_TYPE=RelWithDebInfo

if test -d 'bin'
then
	echo The 'bin' directory exists!
else
	echo Making bin directory . . .
	mkdir bin
fi

echo Building project . . .
cd bin
cmake -DCMAKE_BUILD_TYPE=$BUILD_TYPE .. && make


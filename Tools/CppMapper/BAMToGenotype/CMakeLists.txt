cmake_minimum_required(VERSION 2.8)
project(BAMTOGenotypePrj)

#set(CXX_FLAGS "${CXX_FLAGS} -O3 -Wall -Wextra -pedantic -Wno-unused-function -std=c++20 -m64")
set(CXX_FLAGS "${CXX_FLAGS} -g -Wall -Wextra -pedantic -Wno-unused-function -std=c++20 -m64")

set(LIBS_SRC
	${CMAKE_SOURCE_DIR}/../Libs/logger.cpp
	${CMAKE_SOURCE_DIR}/../Libs/SAM-to-genotypes.cpp)

include_directories(${CMAKE_SOURCE_DIR}/../Libs)

add_library(common_libs ${LIBS_SRC})

add_executable(BAMToGenotype main-bam-to-genotype.cpp)
target_link_libraries(BAMToGenotype common_libs z)

add_executable(UnitTests ../UnitTests/main_unit_tests.cpp)
target_link_libraries(UnitTests common_libs z)


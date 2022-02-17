cmake_minimum_required(VERSION 3.2)
project(FastSTRUCTURE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fPIC")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic")

include_directories(${CMAKE_SOURCE_DIR}/../Libs)

set(VB_SOURCES
	${CMAKE_SOURCE_DIR}/../Libs/logger.cpp
	${CMAKE_SOURCE_DIR}/vb_main.cpp)

add_executable(vb ${VB_SOURCES})


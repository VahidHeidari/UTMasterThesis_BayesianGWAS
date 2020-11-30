#CXX_FLAGS=-Wall -Wextra -pedantic -Wno-unused-function -std=c++20 -m64 -O3 -I../Libs
CXX_FLAGS=-Wall -Wextra -pedantic -Wno-unused-function -std=c++20 -m64 -g -I../Libs

CXX=g++

LIB_OBJS=logger.o SAM-to-genotypes.o
OBJS=$(LIB_OBJS) main-bam-to-genotype.o
EXE_OUT=BAM2Genotype.exe
TEST_OUT=UnitTests.exe



all:
	$(CXX) $(CXX_FLAGS) -c ../Libs/logger.cpp -o logger.o
	$(CXX) $(CXX_FLAGS) -c ../Libs/SAM-to-genotypes.cpp -o SAM-to-genotypes.o
	$(CXX) $(CXX_FLAGS) -c main-bam-to-genotype.cpp -o main-bam-to-genotype.o
	$(CXX) $(CXX_FLAGS) $(OBJS) -o $(EXE_OUT) -lz
	$(CXX) $(CXX_FLAGS) -c ../UnitTests/main_unit_tests.cpp -o main_unit_tests.o
	$(CXX) $(CXX_FLAGS) $(LIB_OBJS) main_unit_tests.o -o $(TEST_OUT) -lz



clean:
	rm -rf *.o $(EXE_OUT) $(TEST_OUT)


CC=g++
CXXFLAGS=-Wall -Wextra -Wpedantic -g -I../source
OBJS=main.o PackedArray.o test-PackedArray.o
SOURCES=main.cpp ../source/PackedArray.cpp test-PackedArray.cpp



all:
	$(CC) $(CXXFLAGS) -c ../source/PackedArray.cpp -o PackedArray.o
	$(CC) $(CXXFLAGS) -c test-PackedArray.cpp -o test-PackedArray.o
	$(CC) $(CXXFLAGS) -c main.cpp -o main.o
	$(CC) $(CXXFLAGS) $(OBJS) -o UnitTest.exe



clean:
	rm -rf *.o UnitTest.exe
	rm -rf *~


CXX_FLAGS=-std=c++2a -Wall -ILibs -Wno-unused -O3
CXX=g++

OBJS=logger.o vb_main.o
OUT_EXE=FastSTRUCTURE.out



all: 
	$(CXX) $(CXX_FLAGS) -c Libs/logger.cpp -o logger.o
	$(CXX) $(CXX_FLAGS) -c FastSTRUCTURE/vb_main.cpp -o vb_main.o
	$(CXX) $(CXX_FLAGS) $(OBJS) -o $(OUT_EXE)



clean:
	rm -f $(OBJS) $(OUT_EXE)


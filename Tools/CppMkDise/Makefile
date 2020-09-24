CXX_FLAGS=-O3 -std=c++2a -Wall -Wno-unused
#CXX_FLAGS=-g -std=c++2a -Wall -Wno-unused



ifdef OS
	STR_OUTNAME=mk_str.exe
	DIS_OUTNAME=mk_dise.exe
else
	STR_OUTNAME=mk_str.out
	DIS_OUTNAME=mk_dise.out
endif



all:
	g++ $(CXX_FLAGS) mk_str.cpp  -o $(STR_OUTNAME)
	g++ $(CXX_FLAGS) mk_dise.cpp -o $(DIS_OUTNAME)
	cp $(STR_OUTNAME) $(DIS_OUTNAME) ../../


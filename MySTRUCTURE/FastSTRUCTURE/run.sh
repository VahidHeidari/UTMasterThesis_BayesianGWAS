#!/bin/bash

IN_FILE=/cygdrive/d/R/fastSTRUCTURE/Datasets/small-3-faststr.txt 
IN_FILE=/cygdrive/d/R/fastSTRUCTURE/Datasets/genos_K2_L100_D50_N100-faststr.txt 
IN_FILE=/cygdrive/d/R/fastSTRUCTURE/Datasets/genos_K2_L1000_D50_N100-faststr.txt 

./build.sh && \
./bin/vb.exe $IN_FILE > log.txt


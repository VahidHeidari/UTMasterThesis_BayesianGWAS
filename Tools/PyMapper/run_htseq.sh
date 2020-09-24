#!/bin/bash

SRR_NAME=SRR1982742
SRR_NAME=SRR1982647


# Finding sub-directory of SRR file.
BASE_PATH=../../GSE68086/Dataset/FastQ
SUB_DIR=HC
if [ ! -f "$BASE_PATH/$SUB_DIR/$SRR_NAME-sorted.sam" ]; then
	SUB_DIR=CRC
	if [ ! -f "$BASE_PATH/$SUB_DIR/$SRR_NAME-sorted.sam" ]; then
		echo 'Could not find input SRR file!'
		exit 1
	fi
fi

SAM_FILE=$BASE_PATH/$SUB_DIR/$SRR_NAME-sorted.sam
COUNT_OUT=$BASE_PATH/$SUB_DIR/$SRR_NAME.count
GTF_FILE=../../hg38-STAR_sample/References/Homo_sapiens.GRCh38.79.gtf


echo START TIME : `date`
htseq-count --mode union --nonunique none --stranded no --type exon --nprocesses 2 --counts_output $COUNT_OUT $SAM_FILE $GTF_FILE
echo 'END TIME   :' `date`


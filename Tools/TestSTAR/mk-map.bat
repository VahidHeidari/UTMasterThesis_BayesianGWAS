@echo off

echo ----------------------------
echo        START mapping
echo ----------------------------

REM add below option for BAM output
REM --outSAMtype        BAM Unsorted				^

call run-star.bat									^
	--runThreadN        1							^
	--genomeDir         SampleGenome				^
	--readFilesIn       SampleReads\reads.fasta		^
	--outFileNamePrefix SampleReads\


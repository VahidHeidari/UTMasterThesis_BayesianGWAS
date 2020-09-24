@echo off

echo -----------------------------------
echo        START making indices
echo -----------------------------------

call run-star.bat										^
	--runThreadN          1								^
	--runMode             genomeGenerate				^
	--genomeDir           SampleGenome					^
	--genomeFastaFiles    SampleGenome\genome.fasta		^
	--sjdbGTFfile         SampleGenome\genome.gtf		^
	--sjdbOverhang        9								^
	--genomeSAindexNbases 3								^
	--outFileNamePrefix SampleGenome\


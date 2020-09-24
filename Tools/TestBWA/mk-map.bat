@echo off

echo -------------------------
echo        BWA mapping
echo -------------------------

call run-bwa.bat mem SampleGenome\genome.fasta SampleReads\reads.fasta > SampleReads\Aligned.out.sam


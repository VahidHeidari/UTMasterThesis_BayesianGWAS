@echo off

del /s *~ *.pyc *.out *.sam *.bam *.tab
del /q SampleGenome\* SampleReads\*
rd SampleGenome SampleReads


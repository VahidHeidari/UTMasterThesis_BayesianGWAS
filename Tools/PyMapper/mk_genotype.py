import os
import sys

import logger
import mk_map
import run_cygwin


IS_DEBUG = False

# Dataset path
#BASE_PATH = 'D:\\Datasets\\GSE68086\\Dataset\\FastQ'
#BASE_PATH = 'G:\\GSE68086\\CRC-sorted'
#BASE_PATH = 'G:\\GSE68086\\HC-sorted'
BASE_PATH = 'H:\\GSE68086\\HC-sorted'
FASTQ_DIR = BASE_PATH
mk_map.TERMINATE_FILE = os.path.join(FASTQ_DIR, 'terminate_file.txt')       # Termination command file

# Tools
SAM_TO_GENO_EXE = 'D:\\Datasets\\Programs\\CppMapper\\x64\\Release\\SAMToGenotype.exe'
SAMTOOLS_EXE = 'D:\\Datasets\\Programs\\samtools-1.10\\samtools.exe'

# References
REF_DIR = 'D:\\Datasets\\hg38-STAR_sample\\References'
ASM_PATH = os.path.join(REF_DIR, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa.bin')



def CollectFileNames():
    file_names = []
    for dirpath, dirnames, filenames in os.walk(FASTQ_DIR):
        for fl in filenames:
            ext = os.path.splitext(fl)[1] 
            if ext != '.sam' or not fl.endswith('-sorted.sam'):
                continue

            file_name = os.path.join(dirpath, fl)
            file_names.append(file_name)
    return file_names



def RunCommand(cmd):
    run_cygwin.RunCommand(['echo'] + cmd if IS_DEBUG else cmd)



def SAMToGenotype(file_name):
    if os.path.isfile(file_name + '.gtp'):
        logger.Log('    Genotype file exists!')
        return

    logger.Log('    SAM to genotype . . .')
    RunCommand([SAM_TO_GENO_EXE, ASM_PATH, file_name, os.path.dirname(file_name)])



def SAMToBAM(file_name):
    if os.path.isfile(file_name + '.bam'):
        logger.Log('    BAM file exits!')
        return

    logger.Log('    Convert SAM to BAM . . .')
    RunCommand([SAMTOOLS_EXE, 'view', '-b', '-o', file_name + '.bam', file_name])



def DeleteTemps(file_name):
    logger.Log('    Delete SAM file . . .')
    if os.path.isfile(file_name):
        os.remove(file_name)



def main():
    logger.Log('\n\n\n========== START ==========')
    logger.Log('FastQ directory -> ' + FASTQ_DIR)
    run_cygwin.AddCygwinPath()

    file_names = CollectFileNames()              # Collect file names.
    for fn in file_names:
        if mk_map.TerminateProcessing():
            logger.Log('---------- Terminated by user! ----------')
            break

        if os.path.isfile(fn + '.gtp') and os.path.isfile(fn + '.bam') and not os.path.isfile(fn):
            logger.Log('The `' + fn + '` has been processed!')
            continue

        logger.Log('*** Processing -> ' + fn)
        SAMToGenotype(fn)
        SAMToBAM(fn)
        DeleteTemps(fn)

    run_cygwin.RemoveCygwinPath()
    logger.Log('\n==========  END  ==========')


if __name__ == '__main__':
    main()


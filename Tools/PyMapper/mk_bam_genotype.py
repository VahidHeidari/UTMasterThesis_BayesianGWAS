import os
import sys

import logger
import mk_map
import run_cygwin



IS_DEBUG = False

# Dataset path
BASE_PATH = 'D:\\Datasets\\GSE68086\\CRC-sorted'
FASTQ_DIR = BASE_PATH

# Tools
BAM_TO_GENO_EXE = 'D:\\Datasets\\Programs\\CppMapper\\BAMToGenotype\\bin\\BAMToGenotype.exe'

# References
REF_DIR = 'D:\\Datasets'
ASM_PATH = os.path.join(REF_DIR, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')



def CollectFileNames():
    file_names = []
    for dirpath, dirnames, filenames in os.walk(FASTQ_DIR):
        for fl in filenames:
            ext = os.path.splitext(fl)[1]
            if ext != '.bam' or not fl.endswith('-sorted.bam'):
                continue

            file_name = os.path.join(dirpath, fl)
            file_names.append(file_name)
    return file_names



def RunCommand(cmd):
    run_cygwin.RunCommand(['echo'] + cmd if IS_DEBUG else cmd)



def BAMToGenotype(file_name):
    if os.path.isfile(file_name + '.gtp'):
        logger.Log('    Genotype file exists!')
        return

    logger.Log('    BAM to genotype . . .')
    RunCommand([BAM_TO_GENO_EXE, ASM_PATH + '.bin', file_name, FASTQ_DIR])



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

        logger.Log('\n\n\n      *** Processing -> ' + fn)
        BAMToGenotype(fn)

    run_cygwin.RemoveCygwinPath()
    logger.Log('\n==========  END  ==========')


if __name__ == '__main__':
    main()


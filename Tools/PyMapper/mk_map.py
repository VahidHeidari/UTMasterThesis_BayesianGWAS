import os
import sys

import logger
import run_cygwin



IS_DEBUG = False

# Dataset path
BASE_PATH = 'D:\\Datasets\\GSE68086\\Dataset\\FastQ'
#FASTQ_DIR = os.path.join(BASE_PATH, 'HC')
FASTQ_DIR = BASE_PATH
TERMINATE_FILE = os.path.join(FASTQ_DIR, 'terminate_file.txt')      # Termination command file

# Tools
GZIP_EXE = 'C:\\PROGRA~1\\Git\\usr\\bin\\gzip.exe'
BWA_EXE = 'D:\\Datasets\\Programs\\bwa-0.7.17\\bwa.exe'
SAMTOOLS_EXE = 'D:\\Datasets\\Programs\\samtools-1.10\\samtools.exe'
FASTQ_DUMP_PATH = 'D:\\Datasets\\GSE68086\\programs\\sratoolkit.2.9.2-win64\\bin\\fastq-dump.exe'
SAM_TO_GENO_EXE = 'D:\\Datasets\\Programs\\CppMapper\\x64\\Release\\SAMToGenotype.exe'

# References
REF_DIR = 'D:\\Datasets\\hg38-STAR_sample\\References'
ASM_PATH = os.path.join(REF_DIR, 'Homo_sapiens.GRCh38.dna.primary_assembly.fa')



def GetBasename(file_name):
    basename = os.path.basename(file_name)
    if os.path.splitext(file_name)[1] == '.gz':
        basename = basename[:-len('.gz')]
    basename = os.path.splitext(basename)[0]
    return basename


def GetBasenamePath(file_name):
    basename = GetBasename(file_name)
    basename_path = os.path.dirname(file_name)
    basename_path = os.path.join(basename_path, basename)
    return basename_path


def MakeFastQPath(file_name):
    return GetBasenamePath(file_name) + '.fastq'

def MakeSAMPath(file_name):
    return GetBasenamePath(file_name) + '.sam'

def MakeSortOutputPath(file_name):
    return GetBasenamePath(file_name) + '-sorted.sam'

def MakeGenotypePath(file_name):
    return GetBasenamePath(file_name) + '-sorted.sam.gtp'

def MakeBAMPath(file_name):
    return GetBasenamePath(file_name) + '-sorted.bam'



def CollectFileNames():
    file_names = []
    for dirpath, dirnames, filenames in os.walk(FASTQ_DIR):
        for fl in filenames:
            ext = os.path.splitext(fl)[1] 
            if ext not in ['.gz', '.sra']:
                continue

            file_name = os.path.join(dirpath, fl)
            if ext == '.gz':
                file_names.insert(0, file_name)
            else:
                file_names.append(file_name)
    return file_names



def RunCommand(cmd):
    run_cygwin.RunCommand(['echo'] + cmd if IS_DEBUG else cmd)



def Decompress(file_name):
    # Make file names.
    out_name = MakeFastQPath(file_name)
    if os.path.isfile(out_name):
        logger.Log('    Decompressed output exists!')
        return

    ext = os.path.splitext(file_name)[1]
    if ext == '.gz':
        logger.Log('    Decompress .tar.gz file . . .')
        RunCommand([GZIP_EXE, '--keep', '--decompress', file_name])
    else:
        logger.Log('    Decompress .sra file . . .')
        RunCommand([FASTQ_DUMP_PATH, file_name, '-O', os.path.dirname(file_name)])
    logger.Log('    Decompression done!')



def Align(file_name):
    # Make file names.
    fastq_name = MakeFastQPath(file_name)
    out_SAM = MakeSAMPath(file_name)
    if os.path.isfile(out_SAM):
        logger.Log('    Align output exists!')
        return

    logger.Log('    Align . . .')
    RunCommand([BWA_EXE, 'mem', '-t', '2', '-o', out_SAM, ASM_PATH, fastq_name])
    logger.Log('    Align done!')



def Sort(file_name):
    # Make file names.
    SAM_path = MakeSAMPath(file_name) 
    sort_path = MakeSortOutputPath(file_name)
    if os.path.isfile(sort_path):
        logger.Log('    Sort output exists!')
        return

    logger.Log('    Sort . . .')
    RunCommand([SAMTOOLS_EXE, 'sort', '-o', sort_path, '-O', 'SAM', SAM_path])
    logger.Log('    Sort done!')



def SAMToGenotype(file_name):
    sorted_sam = MakeSortOutputPath(file_name)
    gtp_path = MakeGenotypePath(file_name)
    if os.path.isfile(gtp_path):
        logger.Log('    Genotype file exists!')
        return

    logger.Log('    SAM to genotype . . .')
    RunCommand([SAM_TO_GENO_EXE, ASM_PATH + '.bin', sorted_sam, os.path.dirname(file_name)])



def SAMToBAM(file_name):
    bam_path = MakeBAMPath(file_name)
    sorted_sam = MakeSortOutputPath(file_name)
    if os.path.isfile(bam_path):
        logger.Log('    BAM file exits!')
        return

    logger.Log('    Convert SAM to BAM . . .')
    RunCommand([SAMTOOLS_EXE, 'view', '-b', '-o', bam_path, sorted_sam])



def DeleteTemps(file_name):
    # Make file names.
    SRR_fastq = MakeFastQPath(file_name)
    SRR_sam = MakeSAMPath(file_name)
    sorted_sam = MakeSortOutputPath(file_name)

    logger.Log('    Delete temps . . .')
    if os.path.isfile(SRR_fastq):
        os.remove(SRR_fastq)
    if os.path.isfile(SRR_sam):
        os.remove(SRR_sam)
    if os.path.isfile(sorted_sam):
        os.remove(sorted_sam)
    logger.Log('    Delete temps done!')



def TerminateProcessing():
    if not os.path.isfile(TERMINATE_FILE):
        return False

    logger.Log('There is a terminate file!')
    with open(TERMINATE_FILE, 'r') as f:
        ln = f.readline().strip()
        logger.Log('file content : `' + ln + '`')
        if ln == 'TERMINATE':
            return True

    logger.Log('There is not any termination command!')
    return False



def main():
    logger.Log('\n\n\n========== START ==========')
    logger.Log('FastQ directory -> ' + FASTQ_DIR)
    run_cygwin.AddCygwinPath()

    file_names = CollectFileNames()              # Collect file names.
    for fn in file_names:
        if TerminateProcessing():
            logger.Log('---------- Terminated by user! ----------')
            break

        if os.path.isfile(MakeGenotypePath(fn)) and os.path.isfile(MakeBAMPath(fn)):
            logger.Log('The `' + fn + '` has been processed!')
            continue

        logger.Log('*** Processing -> ' + fn)
        Decompress(fn)                          # Decompress FastQ gz file.
        Align(fn)                               # Align short reads.
        Sort(fn)                                # Sort SAM file.
        SAMToGenotype(fn)                       # Make genotype file.
        SAMToBAM(fn)                            # Compress SAM file to BAM.
        DeleteTemps(fn)                         # Delete temporary files.

    run_cygwin.RemoveCygwinPath()
    logger.Log('\n==========  END  ==========')


if __name__ == '__main__':
    main()


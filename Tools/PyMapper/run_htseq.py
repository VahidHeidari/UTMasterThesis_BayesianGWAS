import os

import logger
import mk_map
import run_cygwin



IS_DEBUG = False

# Dataset path
LINUX_BASE = '/media/vahidlinux20/0436211A36210E6C/'
BASE_PATH = LINUX_BASE + 'Datasets/GSE68086/Dataset/FastQ'
FASTQ_DIR = BASE_PATH

#LINUX_BASE = '/media/vahidlinux/CA48434048432A91/'
#BASE_PATH = LINUX_BASE + 'GSE68086/CRC-sorted'
#FASTQ_DIR = BASE_PATH

mk_map.TERMINATE_FILE = os.path.join(FASTQ_DIR, 'terminate_file.txt')      # Termination command file

# References
REF_DIR = '/media/vahidlinux20/0436211A36210E6C/' + 'Datasets/hg38-STAR_sample/References'
#REF_DIR = '/media/vahidlinux/0436211A36210E6C/' + 'Datasets/hg38-STAR_sample/References'
GTF_FILE = os.path.join(REF_DIR, 'Homo_sapiens.GRCh38.79.gtf')



def GetBasename(file_name):
    return os.path.basename(file_name)[:-len('-sorted.sam.bam')]

def MakeOutputPath(file_name):
    base_name = GetBasename(file_name)
    count_out = os.path.join(os.path.dirname(file_name), base_name + '.count')
    return count_out



def CollectFileNames():
    file_names = []
    for dirpath, dirnames, filenames in os.walk(FASTQ_DIR):
        for fl in filenames:
            ext = os.path.splitext(fl)[1] 
            if ext != '.bam' or not fl.endswith('-sorted.sam.bam'):
                continue

            file_name = os.path.join(dirpath, fl)
            file_names.append(file_name)
    return file_names



def RunCommand(cmd):
    run_cygwin.RunCommand(['echo'] + cmd if IS_DEBUG else cmd)



def main():
    logger.Log('\n\n\n========== START HTseq ==========')
    logger.Log('FastQ directory -> ' + FASTQ_DIR)
    run_cygwin.AddCygwinPath()

    file_names = CollectFileNames()              # Collect file names.
    for fn in file_names:
        if mk_map.TerminateProcessing():
            logger.Log('---------- Terminated by user! ----------')
            break

        count_out = MakeOutputPath(fn)
        if os.path.isfile(count_out):
            logger.Log('The `' + fn + '` has been processed!')
            continue

        logger.Log('*** Processing -> ' + fn)
        cmd = ['htseq-count', '--mode', 'union', '--nonunique', 'none',
                '--stranded', 'no', '--type', 'exon', '--nprocesses', '4',
                '--counts_output', count_out, fn, GTF_FILE]
        RunCommand(cmd)
        #break

    run_cygwin.RemoveCygwinPath()
    logger.Log('\n==========  END HTseq  ==========')


if __name__ == '__main__':
    main()


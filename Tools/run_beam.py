import os
import subprocess
import sys
import time



IS_DEBUG = False
TEST_NUM = 'test14'
TEST_BASE_PATH = 'F:\\C++\\fastStructure\\TestStr\\Cpp\\' + TEST_NUM + '\\'
BEAM_FILE = 'sm-beam-genos.txt'
BEAM_FILE = 'k1-beam-genos.txt'
#BEAM_FILE = 'k2-beam-genos.txt'

DIS_LOCI = [
    1720, 1721, 1722, 1723,
    1724,
    1725, 1726, 1727, 1728,

    2616, 2617, 2618, 2619,
    2620,
    2621, 2622, 2623, 2624,

    6001, 6002, 6003, 6004,
    6005,
    6006, 6007, 6008, 6009,

    574, 575, 576, 577,
    578,
    579, 580, 581, 582,

    3782, 3783, 3784, 3785,
    3786,
    3787, 3788, 3789, 3790,

    9966, 9967, 9968, 9969,
    9970,
    9971, 9972, 9973, 9974,
]



def Log(msg):
    with open('log.txt', 'a') as f:
        f.write(str(msg) + '\n')
    print(str(msg))



def RunCommand(cmd):
    if IS_DEBUG:
        if os.name == 'nt':
            cmd = ['cmd', '/C', 'echo'] + cmd
        else:
            cmd = ['echo'] + cmd

    cmd = [str(itm) for itm in cmd]
    Log(' START : ' + time.ctime())
    Log(' '.join(cmd))
    ret = subprocess.call(cmd)
    Log('ret : ' + str(ret))
    Log(' END   : ' + time.ctime())
    return ret



def MakeBEAMFile():
    geno_path = TEST_BASE_PATH + 'cpp_dise_genos_K2_L10000_D5_N500.str-faststr.txt'
    label_path = TEST_BASE_PATH + 'cpp_dise_label_K2_L10000_D5_N500.txt'
    beam_path = TEST_BASE_PATH + BEAM_FILE
    if not os.path.isfile(geno_path) or not os.path.isfile(label_path):
        Log('Could not find genotype and label input files!')
        exit(2)

    Log('Genotype input -> ' + geno_path)
    Log('Label input    -> ' + label_path)
    Log('BEAM output    -> ' + beam_path)

    # Read my structure genotype file.
    with open(geno_path, 'r') as f:
        NUM_INDIVS = int(f.readline().split(':')[1].strip())
        NUM_LOCI = int(f.readline().split(':')[1].strip())
        NUM_CLUSTERS = int(f.readline().split(':')[1].strip())
        f.readline()
        TOTAL_INDIVS = NUM_INDIVS * NUM_CLUSTERS
        indivs = [f.readline() for i in range(TOTAL_INDIVS)]

    # Read labels.
    with open(label_path, 'r') as f:
        labels = [1 if f.readline().strip() == 'YES' else 0 for i in range(TOTAL_INDIVS)]

    # Write into BEAM file format.
    with open(beam_path, 'w') as f:
        f.write('ID Chr Pos ')
        #for l in range(len(labels)):
        for l in range(len(labels) // 2):
        #for l in range(len(labels) // 2, TOTAL_INDIVS):
            f.write(str(labels[l]))
            if l + 1 < TOTAL_INDIVS:
                f.write(' ')
        f.write('\n')

        for l in range(NUM_LOCI):
        #for l in DIS_LOCI:
            f.write('rs{} Chr1 {} '.format(l + 1, l))

            #for i in range(TOTAL_INDIVS):
            for i in range(NUM_INDIVS):
            #for i in range(NUM_INDIVS, TOTAL_INDIVS):
                f.write(indivs[i][l])
                if i + 1 < TOTAL_INDIVS:
                    f.write(' ')
            if l + 1 < NUM_LOCI:
                f.write('\n')



def main():
    Log('START : ' + time.ctime())
    indata = TEST_BASE_PATH + BEAM_FILE
    if not os.path.isfile(indata):
        MakeBEAMFile()
    outdata = 'beam1-' + TEST_NUM + '.txt'
    RunCommand(['Release\\BEAMsource.exe', indata, outdata])
    Log('End : ' + time.ctime() + '\n')


if __name__ == '__main__':
    main()


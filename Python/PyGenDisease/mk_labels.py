import os
import random

import numpy as np



# Parameters
NUM_LOCI = 5
NUM_INDIVS = 100
NUM_CHROMOSOMES = 2

# Pathes
FREQS_PATH = 'freqs.txt'
GENOTYPE_PATH = 'genos.txt'
LABELS_FILE = 'labels.txt'



def MultVarSigmoid(params, xs):
    sm = params[0]                  # Bias is the first element in params array.
    for i in range(len(xs)):
        sm += params[i + 1] * xs[i]
    res = 1.0 / (1.0 + np.exp(-sm))
    return res




def MakeRandomFreqs():
    with open(FREQS_PATH, 'w') as f:
        for c in range(NUM_CHROMOSOMES):       # Number of chromosomes is 2.
            for l in range(NUM_LOCI):          # Number of loci
                r = random.uniform(0, 1)
                f.write(str(r))
                if l + 1 != NUM_LOCI:
                    f.write(' ')
            if c + 1 != NUM_CHROMOSOMES:
                f.write('\n')



def GetFreqs():
    if not os.path.isfile(FREQS_PATH):
        print('Making frequencies . . .')
        MakeRandomFreqs()
    with open(FREQS_PATH, 'r') as f:
        freqs = [l.strip().split() for l in f]
    return freqs


def MakeRandomGenotypes():
    freqs = GetFreqs()
    genos = []
    for i in range(NUM_INDIVS):
        geno = []
        for l in range(NUM_LOCI):
            a = 0 if (random.uniform(0, 1) < float(freqs[0][l])) else 1
            b = 0 if (random.uniform(0, 1) < float(freqs[1][l])) else 1
            geno.append(a + b)
        genos.append(geno)

    with open(GENOTYPE_PATH, 'w') as f:
        for i in range(len(genos)):
            geno = genos[i]
            for l in range(NUM_LOCI):
                f.write(str(geno[l]))
                if l + 1 != NUM_LOCI:
                    f.write(' ')
            if i + 1 != NUM_INDIVS:
                f.write('\n')
    return genos



def GetGenotypes():
    if not os.path.isfile(GENOTYPE_PATH):
        print('Makeing genotypes . . .')
        return MakeRandomGenotypes()
    with open(GENOTYPE_PATH, 'r') as f:
        genos = [l.split() for l in f]
    return genos


#
# SNP list format: [ rec, ... , rec ]
#
#   Format of each rec: [ (LOCUS, SNP), ...  ]
#       LOCUS: List of loci from {0, ... , NUM_LOCI - 1}
#       SNP: SNP number is from {0, 1, 2}
#
def GetSNPs(geno, SNPs_loci_lst):
    SNPs = []
    for i in range(len(SNPs_loci_lst)):
        SNPs.append(1)
        rec = SNPs_loci_lst[i]
        for r in rec:
            l = r[0]
            s = r[1]
            SNPs[i] *= (1 if int(geno[l]) == s else 0)
    return SNPs



def MakeLabels(genos):
    if os.path.isfile(LABELS_FILE):
        return

    num_cases = 0
    print('Making labels . . .')
    with open(LABELS_FILE, 'w') as f:
        for g in genos:
            SNPs = GetSNPs(g, [
                [(0, 0)],
                [(2, 1), (3, 1)]
            ])
            p = MultVarSigmoid([-0.1, 0.5, 1], SNPs)
            if p > 0.5:
                num_cases += 1
            lbl = ('YES' if p > 0.5 else 'NO')
            f.write('{}\t{:0.2f}\t{}\n'.format(lbl, p, SNPs))
    print('NUM_CASES: ' + str(num_cases))



def main():
    genos = GetGenotypes()
    lbls = MakeLabels(genos)


if __name__ == '__main__':
    main()


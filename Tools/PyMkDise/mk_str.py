import copy
import random
import time

import mk_my_faststr



NUM_LOCI = 10000
NUM_INDIVS = 100
NUM_CLUSTERS = 2
NUM_CHROMOSOMES = 2

NUM_DIFFS = max(1, int(0.05 * NUM_LOCI + 0.5))

FREQS_PATH = 'freqs_K{}_L{}_D{}.txt'
GENOS_PATH = 'genos_K{}_L{}_D{}_N{}.str'



def MakeAllDiffFreqs():
    freqs = [[[random.uniform(0, 1) for l in range(NUM_LOCI)] for c in range(NUM_CHROMOSOMES)] for k in range(NUM_CLUSTERS)]
    with open('freqs.txt', 'w') as f:
        f.write('NUM_CLUSTERS: {}\n'.format(NUM_CLUSTERS))
        f.write('NUM_LOCI: {}\n'.format(NUM_LOCI))
        f.write('\n')
        for sub in freqs:
            for ch in sub:
                for loc in ch:
                    f.write('{:0.02f} '.format(loc))
                f.write('\n')
            f.write('\n')
    return freqs



def MakeGenos(freqs, num_indivs=NUM_INDIVS):
    genos_path = GENOS_PATH.format(NUM_CLUSTERS, NUM_LOCI, NUM_DIFFS, num_indivs)
    print(time.ctime(), 'Makding -> ' + genos_path + ' . . .')
    with open(genos_path, 'w') as f:
        for k in range(NUM_CLUSTERS):
            for i in range(num_indivs):
                indiv_idx = k * num_indivs + i
                for c in range(NUM_CHROMOSOMES):
                    #f.write('INDIV_{:04d} {} 0 0 0 0\t'.format(indiv_idx + 1, k))       # Structure columns
                    f.write('INDIV_{:04d} '.format(indiv_idx + 1))                      # Structure columns
                    for l in range(NUM_LOCI):
                        r = random.uniform(0, 1)
                        G = 0 if r < freqs[k][c][l] else 1
                        #if indiv_idx == 0:
                        #    print('c', c, 'l', l, 'r', r, 'freq', freqs[k][c][l], 'G', G)
                        f.write(str(G))
                        if l + 1 < NUM_LOCI:
                            f.write(' ')
                    if indiv_idx + 1 < NUM_CLUSTERS * num_indivs or c + 1 < NUM_CHROMOSOMES:
                        f.write('\n')
    mk_my_faststr.MakeMyFastStructure(genos_path, genos_path + '-faststr.txt', NUM_CLUSTERS)



def MakeFreqs():
    print(time.ctime(), 'MakeFreqs()')

    base = [[random.uniform(0, 1) for l in range(NUM_LOCI)] for c in range(NUM_CHROMOSOMES)]
    freqs = [copy.deepcopy(base) for k in range(NUM_CLUSTERS)]

    diffs = []
    for d in range(NUM_DIFFS):
        r = int(random.uniform(0, NUM_LOCI))
        while r in diffs:
            r = int(random.uniform(0, NUM_LOCI))
        diffs.append(r)
    diffs = sorted(diffs)

    for k in range(NUM_CLUSTERS):
        for c in range(NUM_CHROMOSOMES):
            for d in range(NUM_DIFFS):
                freqs[k][c][diffs[d]] = random.uniform(0, 1)
                #if k == 0:
                #    freqs[k][c][diffs[d]] = 0.01 if c == 0 else 0.9
                #else:
                #    freqs[k][c][diffs[d]] = 0.9 if c == 0 else 0.9

    freqs_path = FREQS_PATH.format(NUM_CLUSTERS, NUM_LOCI, NUM_DIFFS)
    with open(freqs_path, 'w') as f:
        f.write('NUM_CLUSTERS: {}\n'.format(NUM_CLUSTERS))
        f.write('NUM_LOCI: {}\n'.format(NUM_LOCI))
        f.write('\n')

        for d in range(NUM_DIFFS):
            f.write('{}'.format(diffs[d]))
            if d + 1 < NUM_DIFFS:
                f.write(' ')
        f.write('\n\n')

        for sub in freqs:
            for ch in sub:
                for loc in ch:
                    f.write('{:0.02f} '.format(loc))
                f.write('\n')
            f.write('\n')
    return freqs



def main():
    print(time.ctime(), 'num_diffs', NUM_DIFFS, 'of', NUM_LOCI)
    freqs = MakeFreqs()
    #genos = MakeGenos(freqs, 10000)
    genos = MakeGenos(freqs, 1000)
    #genos = MakeGenos(freqs, 500)
    #genos = MakeGenos(freqs, 100)



if __name__ == '__main__':
    main()


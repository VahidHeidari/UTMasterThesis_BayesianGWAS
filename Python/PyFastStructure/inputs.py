import os
import random
import struct


DATA_SET = 'genos.txt'
#BINARY_DATA_SET = 'F:\\Python\\PyGenDisease\\Genos\\genos_K2_L100_D0.5.bin'
BINARY_DATA_SET = 'F:\\Python\\PyGenDisease\\Genos\\genos_K2_L200_D0.02.bin'

NUM_CLUSTERS = 2
NUM_INDIVS = 100
NUM_CHROMOSOMES = 2
NUM_LOCI = 100


def MakeRandomAlleleFrequencies():
    freqs = [[[random.uniform(0, 1) for l in xrange(NUM_LOCI)] for c in xrange(NUM_CHROMOSOMES)] for k in xrange(NUM_CLUSTERS)]
    with open('freqs.txt', 'w') as f:
        f.write('NUM_CLUSTERS:{}\n'.format(NUM_CLUSTERS))
        f.write('NUM_INDIVS:{}\n'.format(NUM_INDIVS))
        f.write('NUM_CHROMOSOMES:{}\n'.format(NUM_CHROMOSOMES))
        f.write('NUM_LOCI:{}\n'.format(NUM_LOCI))
        f.write('\n')
        for clusters in freqs:
            for chromosomes in clusters:
                for locus in chromosomes:
                    f.write('{:.2f} '.format(locus))
                f.write('\n')
            f.write('\n')
    return freqs


def MakeRandomGenotypesFile(freqs, num_indivs):
    K = len(freqs)
    C = len(freqs[0])
    L = len(freqs[0][0])
    with open(DATA_SET, 'w') as f:
        for k in xrange(K):
            freq = freqs[k]
            for ind in xrange(num_indivs):
                for chrom in xrange(C):
                    c = freq[chrom]
                    for loc in xrange(L):
                        allele_freq = c[loc]
                        allele = 0 if random.uniform(0, 1) < allele_freq else 1
                        f.write('{} '.format(allele))
                    f.write('\n')
                f.write('\n')
            f.write('\n\n')


def MakeRandomGenotypes():
    print('MakeRandomGenotypes')
    freqs = MakeRandomAlleleFrequencies()
    MakeRandomGenotypesFile(freqs, NUM_INDIVS)


def ReadGenotypes():
    if not os.path.isfile(DATA_SET):
        MakeRandomGenotypes()

    indivs = []
    chrom = 0
    ind = []
    with open(DATA_SET, 'r') as f:
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue

            chromosome = [int(s) for s in line.split(' ')]
            ind.append(chromosome)
            chrom += 1
            if chrom == NUM_CHROMOSOMES:
                indivs.append(ind)
                chrom = 0
                ind = []
    return indivs


def ReadBinaryGenotypes():
    with open(BINARY_DATA_SET, 'rb') as f:
        hdr = struct.unpack('<III', f.read(3 * 4))
        num_clusters = hdr[0]
        num_loci = hdr[1]
        num_indivs = hdr[2]
        indivs = []
        for i in range(num_indivs):
            ind = []
            bit_num = 0
            buff = struct.unpack('<B', f.read(1))[0]
            for l in range(num_loci):
                if bit_num == 8:
                    buff = struct.unpack('<B', f.read(1))[0]
                    bit_num = 0
                g = (buff >> bit_num) & 0x03
                ind.append(g)
                bit_num += 2
            indivs.append(ind)

    data_set = []
    for i in range(num_indivs):
        ind = [[], []]
        for l in range(num_loci):
            g = indivs[i][l]
            if g == 0:
                ind[0].append(0)
                ind[1].append(0)
            elif g == 1:
                ind[0].append(1)
                ind[1].append(0)
            else:
                ind[0].append(1)
                ind[1].append(1)
        data_set.append(ind)
    return data_set


def GetNumIndivs(data_set):
    return len(data_set)


def GetNumChr(data_set):
    return len(data_set[0])


def GetNumLoci(data_set):
    return len(data_set[0][0])


def GetIndiv(data_set, indiv_idx):
    return data_set[indiv_idx]


def GetChromosome(data_set, indiv_idx, chromosome_idx):
    return GetIndiv(data_set, indiv_idx)[chromosome_idx]


def GetAllele(data_set, indiv_idx, chromosome_idx, locus_idx):
    return GetChromosome(data_set, indiv_idx, chromosome_idx)[locus_idx]


# ---------- G ----------
def MakeGenotypes(data_set):
    G = []
    for i in xrange(GetNumIndivs(data_set)):
        gene = [0 for l in xrange(GetNumLoci(data_set))]
        for l in xrange(GetNumLoci(data_set)):
            for c in xrange(GetNumChr(data_set)):
                gene[l] += GetAllele(data_set, i, c, l)
        G.append(gene)
    return G


def GetGenotype(genotypes, individual_idx, locus_idx):
    return genotypes[individual_idx][locus_idx]


if __name__ == '__main__':
    print('This is a library :)')

import itertools
import os
import random
import struct
import time

import conf_file
import mk_pop_freqs



BASE_FREQS_PATH = 'PopConfigs\\base_freqs.bin'
CONF_PATH_FMT = 'PopConfigs\\sub_freqs_L{}.bin'
GENO_PATH_FMT = 'Genos\\genos_N{}_K{}_L{}_D{}.bin'
NUM_INDIVS = 1000
MAX_NUM_LOCI = 100000
NUM_CHROMOSOMES = 2
DIFF_MAF_LIST = [0.5] #[0.02] #[0.02, 0.05, 0.10, 0.20]
NUM_POP_LIST = [2] #[2, 3, 4, 5]
LOCI_LIST = [100] #[i * 100 for i in range(1, 10)] + [i * 1000 for i in range(1, 101)]



def FloatToByte(f):
    return int(f * 100)

def ByteToFloat(b):
    return float(b) / 100.0



def MakeBaseFreqs():
    print('Making Base Freqs . . .')
    with open(BASE_FREQS_PATH, 'w') as f:
        for c in xrange(NUM_CHROMOSOMES):
            for l in xrange(MAX_NUM_LOCI):
                f.write('{:0.03f}'.format(random.uniform(0, 1)))
                if l + 1 < MAX_NUM_LOCI:
                    f.write(' ')
            if c + 1 < NUM_CHROMOSOMES:
                f.write('\n')



def MakeBinaryBaseFreqs():
    print('Making Base Freqs . . .')
    with open(BASE_FREQS_PATH, 'wb') as f:
        f.write(struct.pack('<I', MAX_NUM_LOCI))
        for c in xrange(NUM_CHROMOSOMES):
            for l in xrange(MAX_NUM_LOCI):
                f.write(struct.pack('<B', FloatToByte(random.uniform(0, 1))))



def ReadBinaryBaseFreqs():
    base_freqs = [[] for c in xrange(NUM_CHROMOSOMES)]
    with open(BASE_FREQS_PATH, 'rb') as f:
        num_loci = struct.unpack('<I', f.read(4))[0]
        for c in xrange(NUM_CHROMOSOMES):
            for l in xrange(num_loci):
                p = ByteToFloat(struct.unpack('<B', f.read(1))[0])
                base_freqs[c].append(p)
    return base_freqs



def MakeDiffMAFs(num_loci):
    old_diffs = []
    for d in DIFF_MAF_LIST:
        conf = conf_file.ConfFile()
        conf.key_vals = {'NUM_LOCI' : num_loci, 'DIFF_MAFS' : d}
        all_diffs = list(itertools.chain.from_iterable(old_diffs))
        diffs = mk_pop_freqs.SelectDiffMAFs(conf, all_diffs)
        old_diffs.append(diffs)
    return old_diffs



def MakeSubFreqs(num_loci):
    num_diffs = max(1, int(max(DIFF_MAF_LIST) * num_loci + 0.5))
    num_pops = max(NUM_POP_LIST)
    sub_freqs = []
    for k in range(num_pops):
        freqs = []
        for c in range(NUM_CHROMOSOMES):
            freqs.append([random.uniform(0, 1) for p in range(num_diffs)])
        sub_freqs.append(freqs)
    return sub_freqs



def DumpConfs(conf_path, diffs, sub_freqs):
    NUM_POPS = len(sub_freqs)
    NUM_SUB_FREQS_LOCI = len(sub_freqs[0][0])
    with open(conf_path, 'w') as f:
        # MAFs
        for i in range(len(diffs)):
            rec = diffs[i]
            for d in range(len(rec)):
                f.write(str(diffs[i][d]))
                if d + 1 < len(rec):
                    f.write(' ')
            f.write('\n')
        f.write('\n\n')

        # Sub-population frequenceis.
        for k in range(NUM_POPS):
            for c in range(NUM_CHROMOSOMES):
                for l in range(NUM_SUB_FREQS_LOCI):
                    f.write('{:0.03f}'.format(sub_freqs[k][c][l]))
                    if l + 1 < NUM_SUB_FREQS_LOCI:
                        f.write(' ')
                f.write('\n')
            if k + 1 < NUM_POPS:
                f.write('\n')



def DumpBinaryConfs(conf_path, num_loci, diffs, sub_freqs):
    NUM_POPS = len(sub_freqs)
    NUM_DIFFS = len(diffs)
    NUM_SUB_FREQS_LOCI = len(sub_freqs[0][0])
    with open(conf_path, 'wb') as f:
        # Header
        f.write(struct.pack('<BBII', NUM_POPS, NUM_DIFFS, NUM_SUB_FREQS_LOCI, num_loci))
        for d in range(len(diffs)):
            f.write(struct.pack('<I', len(diffs[d])))

        # MAFs
        for i in range(len(diffs)):
            for d in range(len(diffs[i])):
                if num_loci <= 65000:
                    f.write(struct.pack('<H', diffs[i][d]))
                else:
                    f.write(struct.pack('<I', diffs[i][d]))

        # Sub-population frequenceis.
        for k in range(NUM_POPS):
            for c in range(NUM_CHROMOSOMES):
                for l in range(NUM_SUB_FREQS_LOCI):
                    f.write(struct.pack('<B', FloatToByte(sub_freqs[k][c][l])))



def ReadBinaryConf(conf_path):
    with open(conf_path, 'rb') as f:
        # Header
        hdr = struct.unpack('<BBII', f.read(10))
        NUM_POPS = hdr[0]
        NUM_DIFFS = hdr[1]
        NUM_SUB_FREQS_LOCI = hdr[2]
        num_loci = hdr[3]

        diffs_len = []
        for d in range(NUM_DIFFS):
            dln = struct.unpack('<I', f.read(4))[0]
            diffs_len.append(dln)

        # MAFs
        diffs = [[] for i in range(NUM_DIFFS)]
        for i in range(NUM_DIFFS):
            for d in range(diffs_len[i]):
                if num_loci <= 65000:
                    diffs[i].append(struct.unpack('<H', f.read(2))[0])
                else:
                    diffs[i].append(struct.unpack('<I', f.read(4))[0])

        # Sub-population frequenceis.
        sub_pops = [ [ [] for c in range(NUM_CHROMOSOMES) ] for j in range(NUM_POPS) ]
        for k in range(NUM_POPS):
            for c in range(NUM_CHROMOSOMES):
                for l in range(NUM_SUB_FREQS_LOCI):
                    p = ByteToFloat(struct.unpack('<B', f.read(1))[0])
                    sub_pops[k][c].append(p)
    return diffs, sub_pops



def GenerateIndividuals(num_pops, num_loci, base_freqs, sub_freqs, diffs):
    indivs = []
    for k in range(num_pops):
        for i in range(NUM_INDIVS):
            geno = []
            d_idx = 0
            for l in range(num_loci):
                if d_idx < len(diffs) and l == diffs[d_idx]:
                    p0 = sub_freqs[k][0][d_idx]
                    p1 = sub_freqs[k][1][d_idx]
                    d_idx += 1
                else:
                    p0 = base_freqs[0][l]
                    p1 = base_freqs[1][l]

                a = 0 if (random.uniform(0, 1) < p0) else 1
                b = 0 if (random.uniform(0, 1) < p1) else 1
                geno.append(a + b)
            indivs.append(geno)
    return indivs



def DumpIndividuals(geno_path, num_pops, num_loci, indivs):
    with open(geno_path, 'wb') as f:
        f.write(struct.pack('<III', num_pops, num_loci, len(indivs)))
        buff = bit_pos = 0
        for i in range(len(indivs)):
            for l in range(num_loci):
                if bit_pos == 8:
                    f.write(struct.pack('<B', buff))
                    buff = bit_pos = 0
                buff |= indivs[i][l] << bit_pos
                bit_pos += 2
            if bit_pos:
                f.write(struct.pack('<B', buff))



def GenerateAndDumpBinaryIndivduals(num_loci, base_freqs):
    for k in NUM_POP_LIST:
        for d in DIFF_MAF_LIST:
            geno_path = GENO_PATH_FMT.format(NUM_INDIVS, k, num_loci, d)
            if os.path.isfile(geno_path):
                continue

            print('Making genotypes {} . . .'.format(geno_path))
            conf_path = CONF_PATH_FMT.format(num_loci)
            diffs, sub_freqs = ReadBinaryConf(conf_path)
            for d in range(len(diffs)):
                all_diffs = list(itertools.chain.from_iterable(diffs[0 : d + 1]))
                indivs = GenerateIndividuals(k, num_loci, base_freqs, sub_freqs, all_diffs)
                DumpIndividuals(geno_path, k, num_loci, indivs)



def main():
    print('START TIME : {}'.format(time.ctime()))
    if not os.path.isfile(BASE_FREQS_PATH):
        MakeBinaryBaseFreqs()
    base_freqs = ReadBinaryBaseFreqs()

    for l in LOCI_LIST:
        conf_path = CONF_PATH_FMT.format(l)
        if os.path.isfile(conf_path):
            continue

        print('Making config for L{} . . .'.format(l))
        diffs = MakeDiffMAFs(l)
        sub_freqs = MakeSubFreqs(l)
        DumpBinaryConfs(conf_path, l, diffs, sub_freqs)

    for l in LOCI_LIST:
        GenerateAndDumpBinaryIndivduals(l, base_freqs)
    print('END   TIME : {}'.format(time.ctime()))



if __name__ == '__main__':
    main()


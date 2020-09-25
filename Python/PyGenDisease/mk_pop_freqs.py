import copy
import struct
import sys
import random

import conf_file



BINARY_HEADER_FMT = '<BBIII'



def MakeRandomFrequencies(conf):
    base_freqs = []
    for c in range(int(conf.key_vals['NUM_CHROMOSOMES'])):
        freq = []
        for l in range(int(conf.key_vals['NUM_LOCI'])):
            freq.append(random.uniform(0, 1))
        base_freqs.append(freq)
    return base_freqs



def SelectDiffMAFs(conf, old_diffs=[]):
    num_diff_MAFs = max(1, int(float(conf.key_vals['DIFF_MAFS']) * float(conf.key_vals['NUM_LOCI']) + 0.5))
    num_diff_MAFs -= len(old_diffs)
    if num_diff_MAFs <= 0:
        return []

    diff_MAFs = []
    for d in range(num_diff_MAFs):
        l = int(random.uniform(0, int(conf.key_vals['NUM_LOCI'])))
        while l in diff_MAFs or l in old_diffs:
            l = int(random.uniform(0, int(conf.key_vals['NUM_LOCI'])))
        diff_MAFs.append(l)
    return sorted(diff_MAFs)



def MakeSubPopulationMAFs(conf, base_freqs, diff_MAFs):
    sub_freqs = []
    NUM_CHROMOSOMES = int(conf.key_vals['NUM_CHROMOSOMES'])
    for k in range(int(conf.key_vals['NUM_POPS'])):
        sub_freqs.append([[0.0 for i in range(len(diff_MAFs))] for c in range(NUM_CHROMOSOMES)])
        for c in range(NUM_CHROMOSOMES):
            for d in range(len(diff_MAFs)):
                sub_freqs[k][c][d] = random.uniform(0, 1)
    return sub_freqs



def DumpFrequencies(conf, diff_MAFs, base_freqs, sub_freqs):
    print('Dump to -> ' + conf.key_vals['OUTPUT_FREQS'])
    with open(conf.key_vals['OUTPUT_FREQS'], 'w') as f:
        # Configs.
        f.write('# NUM_POPS:\t\t{}\n'.format(conf.key_vals['NUM_POPS']))
        f.write('# NUM_LOCI:\t\t{}\n'.format(conf.key_vals['NUM_LOCI']))
        f.write('# NUM_INDIVS:\t{}\n'.format(conf.key_vals['NUM_INDIVS']))
        f.write('# DIFF_MAFS:\t{}\n'.format(conf.key_vals['DIFF_MAFS']))
        f.write('# OUTPUT_FREQS:\t{}\n'.format(conf.key_vals['OUTPUT_FREQS']))
        f.write('\n\n')

        # Dump base frequenceis.
        f.write('# Base frequencies:\n')
        for c in range(int(conf.key_vals['NUM_CHROMOSOMES'])):
            for l in range(int(conf.key_vals['NUM_LOCI'])):
                f.write('{:0.03f}'.format(base_freqs[c][l]))
                if l + 1 < int(conf.key_vals['NUM_LOCI']):
                    f.write(' ')
            f.write('\n')
        f.write('\n\n')

        # Dump different MAFs.
        f.write('# Diff MAFs:\n')
        for d in diff_MAFs:
            f.write(str(d))
            if d != diff_MAFs[len(diff_MAFs) - 1]:
                f.write(' ')
        f.write('\n\n\n')

        # Dump sub-population frequencies.
        f.write('# Sub-population frequencies:\n')
        for k in range(int(conf.key_vals['NUM_POPS'])):
            for c in range(int(conf.key_vals['NUM_CHROMOSOMES'])):
                for l in range(len(diff_MAFs)):
                    f.write('{:0.03f}'.format(sub_freqs[k][c][l]))
                    if l + 1 < len(diff_MAFs):
                        f.write(' ')
                f.write('\n')
            if k + 1 < int(conf.key_vals['NUM_POPS']):
                f.write('\n')



def DumpBinaryFrequencies(conf, diff_MAFs, base_freqs, sub_freqs):
    f = open(conf.key_vals['OUTPUT_FREQS'], 'wb')

    # Write configs.
    hdr = struct.pack(BINARY_HEADER_FMT,
        conf.key_vals['NUM_POPS'],
        conf.key_vals['NUM_CHROMOSOMES'],
        conf.key_vals['NUM_LOCI'],
        len(diff_MAFs), # Instead of float 'conf.key_vals['DIFF_MAFS']' write integer
        conf.key_vals['NUM_INDIVS'])
    f.write(hdr)

    # Write different MAFs.
    for c in range(int(conf.key_vals['NUM_CHROMOSOMES'])):
        for d in diff_MAFs:
            f.write(struct.pack('<I', d))

    # Write base frequencies.
    for c in range(len(base_freqs)):
        for p in base_freqs[c]:
            f.write(struct.pack('<f', p))

    # Write sub-populations frequencies.
    for k in range(len(sub_freqs)):
        for c in range(len(sub_freqs[0])):
            for p in sub_freqs[k][c]:
                f.write(struct.pack('<f', p))
    f.close()



def MakeSubPopulationConfig(conf):
    base_freqs = MakeRandomFrequencies(conf)                        # Make population allele frequencies.
    diff_MAFs = SelectDiffMAFs(conf)                                # Select loci for diff MAFs sub-populations.
    sub_freqs = MakeSubPopulationMAFs(conf, base_freqs, diff_MAFs)  # Make sub-population MAFs
    DumpFrequencies(conf, diff_MAFs, base_freqs, sub_freqs)         # Dump frequencies.



def main():
    conf_path = sys.argv[1]
    conf = conf_file.ConfFile()
    conf.InitFromFile(conf_path)
    conf.Print()
    MakeSubPopulationConfig(conf)



if __name__ == '__main__':
    main()


import random
import struct

import conf_file
import mk_pop_freqs



BINARY_HEADER_FMT = '<BBIII'



def GetSubFreqs(conf):
    base_freqs = mk_pop_freqs.MakeRandomFrequencies(conf)                        # Make population allele frequencies.
    diff_MAFs = mk_pop_freqs.SelectDiffMAFs(conf)                                # Select loci for diff MAFs sub-populations.
    sub_freqs = mk_pop_freqs.MakeSubPopulationMAFs(conf, base_freqs, diff_MAFs)  # Make sub-population MAFs
    return diff_MAFs, base_freqs, sub_freqs



def MakeConfig(num_pops, num_loci, diff_MAFs):
    # Make configurations.
    out_file = 'PopConfigs\\freqs_K{}_L{}_D{}.bin'.format(num_pops, num_loci, diff_MAFs)
    conf = conf_file.ConfFile()
    conf.key_vals = {
        'NUM_POPS': 	    num_pops,
        'NUM_LOCI': 	    num_loci,
        'DIFF_MAFS':	    diff_MAFs,
        'OUTPUT_FREQS':     out_file,
        'NUM_INDIVS':	    100,
        'NUM_CHROMOSOMES':  2,
    }
    #conf.Print()
    return conf



def MakeIndividuals(conf, sub_freqs):
    genotypes = []
    for k in xrange(int(conf.key_vals['NUM_POPS'])):
        for i in xrange(int(conf.key_vals['NUM_INDIVS'])):
            genotype = []
            for l in xrange(int(conf.key_vals['NUM_LOCI'])):
                a = 0 if (random.uniform(0, 1) < float(sub_freqs[k][0][l])) else 1
                b = 0 if (random.uniform(0, 1) < float(sub_freqs[k][1][l])) else 1
                genotype.append(a + b)
            genotypes.append(genotype)
    return genotypes


def ReadBinaryIndividuals(bin_path):
    f = open(bin_path, 'rb')
    blocks = f.read()

    # Unpack header and init confs.
    hdr_size = struct.calcsize(BINARY_HEADER_FMT)
    fields = struct.unpack(BINARY_HEADER_FMT, blocks[0:hdr_size])
    conf = conf_file.ConfFile()
    conf.key_vals = {
        'NUM_POPS': 	    fields[0],
        'NUM_CHROMOSOMES':  fields[1],
        'NUM_LOCI': 	    fields[2],
        'DIFF_MAFS':	    fields[3],
        'NUM_INDIVS':	    fields[4],
        'OUTPUT_FREQS':     bin_path,
    }

    # Unpack MAFs.
    diff_MAFs = []
    for d in range(conf.key_vals['DIFF_MAFS']):
        blks = blocks[hdr_size + d : hdr_size + d + 4]
        i = struct.unpack('<I', blks)
        diff_MAFs.append(i[0])

    #for c in range(conf.key_vals['NUM_CHROMOSOMES']):
    f.close()
    return conf, diff_MAFs



def DumpIndivs(conf, diff_MAFs, base_freqs, sub_freqs, indivs):
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



def MakeAndDumpIndividuals(num_pops, num_loci, diff_MAFs):
    conf = MakeConfig(num_pops, num_loci, diff_MAFs)
    print('Dump to -> ' + conf.key_vals['OUTPUT_FREQS'])
    diff_MAFs, base_freqs, sub_freqs = GetSubFreqs(conf)
    mk_pop_freqs.DumpBinaryFrequencies(conf, diff_MAFs, base_freqs, sub_freqs)

    #indivs = MakeIndividuals(conf, sub_freqs)
    #DumpIndivs(conf, diff_MAFs, base_freqs, sub_freqs, indivs)

    # Test writing.
    #conf2, diff_MAFs = ReadBinaryIndividuals(conf.key_vals['OUTPUT_FREQS'])
    #conf2.Print()



def main():
    LOCI_LIST = [i * 100 for i in range(1, 10)] + [i * 1000 for i in range(1, 101)]
    for k in [2, 3, 4, 5]:
        for d in [0.02]: #[0.02, 0.05, 0.10, 0.20]:
            for l in LOCI_LIST:
                MakeAndDumpIndividuals(k, l, d)



if __name__ == '__main__':
    main()


import math
import scipy.stats as st

import mk_str



PSEUDOCOUNT = 1e-7



def ReadFreqs():
    all_freqs = []
    with open(mk_str.FREQS_PATH, 'r') as f:
        num_clusters = int(f.readline().split(':')[1].strip())
        num_loci = int(f.readline().split(':')[1].strip())

        f.readline()
        diffs = [int(d) for d in f.readline().split()]

        for l in f:
            l = l.strip()
            if len(l) == 0:
                continue
            all_freqs.append([float(p) for p in l.split()])

    freqs = []
    for i in range(len(all_freqs) / 2):
        freqs.append([all_freqs[i * 2], all_freqs[i * 2 + 1]])
    return num_clusters, num_loci, diffs, freqs



def ReadGenos():
    all_genos = []
    with open(mk_str.GENOS_PATH, 'r') as f:
        for l in f:
            all_genos.append([int(g) for g in l.split()[6:]])

    genos = []
    for i in range(len(all_genos) / 2):
        indiv = []
        for l in range(len(all_genos[0])):
            G = all_genos[i * 2][l] + all_genos[i * 2 + 1][l]
            indiv.append(G)
        genos.append(indiv)
    return genos



def CountGenos(genos, locus, is_normalized_output=True):
    genos_count = [ [0, 0, 0], [0, 0, 0] ]
    for i in range(len(genos) / 2):
        G = genos[i][locus]
        genos_count[0][G] += 1

    for i in range(len(genos) / 2, len(genos)):
        G = genos[i][locus]
        genos_count[1][G] += 1

    if is_normalized_output:
        for i in range(3):
            genos_count[0][i] /= len(genos) / 2.0
            genos_count[1][i] /= len(genos) / 2.0
    return genos_count



def PrintDiffInfos(freqs, locus, cluster_idx):
    pa = freqs[cluster_idx][0][locus]
    pb = freqs[cluster_idx][1][locus]
    G_0 = pa * pb
    G_1 = pa * (1.0 - pb) + pb * (1.0 - pa)
    G_2 = (1.0 - pa) * (1.0 - pb)
    print('freqs[{}][{}]:'.format(cluster_idx, locus))
    print('pa', pa, 'pb', pb)
    print('G_0', G_0, 'G_1', G_1, 'G_2', G_2)
    print('')



def main():
    num_clusters, num_loci, diffs, freqs = ReadFreqs()
    genos = ReadGenos()
    print('diffs', diffs)
    print('infos of diff[0]', diffs[0])
    print('')

    PrintDiffInfos(freqs, diffs[0], 0)
    PrintDiffInfos(freqs, diffs[0], 1)

    df = 0
    sm = 0.0
    for l in range(num_loci):
        genos_count = CountGenos(genos, l)
        #print('locus : ', l)
        #print('Exp:', genos_count[0])
        #print('Obs:', genos_count[1])

        for i in range(len(genos_count[0])):
            E = genos_count[0][i]
            O = genos_count[1][i]
            if E == 0:
                continue

            sm += math.pow(O - E, 2) / E
            df += 1

    p_val = 1 - st.chi2.cdf(sm, df - 1)
    print('Chi2', sm, 'df', df, 'p-value', p_val)
    print('----')


if __name__ == '__main__':
    main()


import time



def MakeMyFastStructure(in_path, out_path, num_clusters):
    print(time.ctime(), 'Makeing MyFastStructure -> ', out_path)

    indivs = []
    with open(in_path, 'r') as f:
        for l in f:
            indivs.append([int(g) for g in l.split()[1:]])

    NUM_INDIVS = len(indivs) / 2
    NUM_LOCI = len(indivs[0])
    print('NUM_INDIVS', NUM_INDIVS, 'NUM_LOCI', NUM_LOCI)

    with open(out_path, 'w') as f:
        f.write('NUM_INDIVS: {}\n'.format(NUM_INDIVS))
        f.write('NUM_LOCI: {}\n'.format(NUM_LOCI))
        f.write('NUM_CLUSTERS: {}\n'.format(num_clusters))
        f.write('\n')
        for i in range(NUM_INDIVS):
            for l in range(NUM_LOCI):
                G = indivs[i * 2][l] + indivs[i * 2 + 1][l]
                f.write(str(G))
            f.write('\n')


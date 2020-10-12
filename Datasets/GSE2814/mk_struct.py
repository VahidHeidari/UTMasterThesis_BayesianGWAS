import read_responses



SNP_INDIC_FILE = 'SNP-indicators.txt'



def IsDisease(weight):
    weight = float(weight)
    #return weight < 4
    #return weight > 6
    #return weight > 5
    #return weight > 4.5
    #return weight > 4.4
    #return weight > 4.3
    return weight > 4.2
    #return weight < 3.9 or weight > 4.3
    #return weight < 4 or weight > 4.2



def ReadSNPIndicators():
    snps = [ [sp.strip() for sp in l.strip().split()] for l in open(SNP_INDIC_FILE, 'r') ]
    snps[0].insert(0, '')
    for i in range(len(snps[0])):
        snps[0][i] = snps[0][i].replace('"', '')
    for i in range(1, len(snps)):
        snps[i][0] = 'F2_' + snps[i][0].replace('"', '')
    return snps



def SNPsToStructure(snps):
    indivs = []
    for i in range(1, len(snps)):
        j = 1
        mice_id = snps[i][0]
        genos = [ [ mice_id ], [ mice_id ] ]
        while j < len(snps[i]):
            g = snps[i][j]
            c = 0 if j % 2 else 1
            genos[c].append(g)
            j += 1
        indivs.append(genos)
    return indivs



def DumpStructure(structure):
    NUM_INDIVS = len(structure)
    NUM_LOCI = len(structure[0][0])
    with open('mice-str.txt', 'wb') as f:
        for i in range(NUM_INDIVS):
            for c in range(2):
                for l in range(NUM_LOCI):
                    f.write(structure[i][c][l])
                    if l + 1 < NUM_LOCI:
                        f.write(' ')
                f.write('\n')



def DumpBEAM(structure, flt_resps):
    NUM_INDIVS = len(structure)
    NUM_LOCI = len(structure[0][0])
    print(len(structure), len(flt_resps))

    f = open('BEAM-genos.txt', 'w')
    f.write('ID Chr Position ')
    empty_idxs = []
    for l in range(len(flt_resps)):
        mid = flt_resps[l][0]
        weight = flt_resps[l][read_responses.RESP_RECS['TotalFat']]
        if weight == '':
            print('empty!', l, mid)
            empty_idxs.append(l)
            continue

        weight = float(weight)
        f.write('0' if IsDisease(weight) else '1')
        if l + 1 < len(flt_resps):
            f.write(' ')
    f.write('\n')
    print(len(empty_idxs))

    for l in range(1, NUM_LOCI):
        f.write('rs{} ch1 {} '.format(l, l))
        for i in range(NUM_INDIVS):
            if i in empty_idxs:
                continue

            ind = structure[i]
            G = [ '00', '01', '10' ].index(ind[0][l] + ind[1][l])
            f.write(str(G))
            if i + 1 < NUM_INDIVS:
                f.write(' ')
        if l + 1 < NUM_LOCI:
            f.write('\n')
    f.close()



def FindIdx(mid, resps):
    for i in range(1, len(resps)):
        if resps[i][read_responses.RESP_RECS['MiceID']] == mid:
            return i
    return None


def FilterResponses(resps, snps):
    fltr = []
    for i in range(1, len(snps)):
        idx = FindIdx(snps[i][0], resps)
        if not idx:
            print('Not found:', i, snps[i][0])
            fltr.append([])
        else:
            fltr.append(resps[idx])
    return fltr



def MakeLabels(flt_resps):
    with open('filtered-Labels.txt', 'w') as f:
        for i in range(len(flt_resps)):
            rec = flt_resps[i]
            weight = rec[8]
            if len(weight):
                lbl = 'NO ' if IsDisease(weight) else 'YES'
            else:
                lbl = '---'
            f.write('{}    {}  {}\n'.format(lbl, rec[0], rec[8]))



def main():
    # Read SNPs.
    snps = ReadSNPIndicators()
    print('Num SNPs rec : ', len(snps))
    print('Num loci x 2 : ', len(snps[1][1:]))
    print('Num loci     : ', len(snps[1][1:]) / 2)

    # Convert to structure format.
    structure = SNPsToStructure(snps)
    print('')
    print('Num Indivs : ', len(structure))
    print('Num Chroms : ', len(structure[0]))
    print('Num loci   : ', len(structure[0][0]))

    # Dump structure formatted file.
    #DumpStructure(structure)

    # Read responses and make labels.
    resps = read_responses.ReadResponses()
    flt_resps = FilterResponses(resps, snps)
    MakeLabels(flt_resps)
    DumpBEAM(structure, flt_resps)
    print('')
    print('Num resps      : ', len(resps))
    print('Filtered resps : ', len(flt_resps))


if __name__ == '__main__':
    main()


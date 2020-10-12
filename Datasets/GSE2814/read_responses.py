import random

import numpy as np
import matplotlib.pyplot as plt



THERESHOLD = 4.3

RESPONSE_FILE = 'Mice and responses.csv'
SAMPLE_SOURCE_FILE = 'Mice and Genes.csv'

SOURCE_NAMES_IDX = 68
ID_REF_IDX = 69



RESP_RECS = {
    'MiceID' : 0,
    'Strain' : 1,
    'sex' : 2,
    'parents' : 3,
    'WeightG' : 4,
    'LengthCM' : 5,
    'AbFat' : 6,
    'OtherFat' : 7,
    'TotalFat' : 8,
    '100xfat/weight' : 9,
    'Trigly' : 10,
    'TotalChol' : 11,
    'HDLChol' : 12,
    'UC' : 13,
    'FFA' : 14,
    'Glucose' : 15,
    'LDL+VLDL' : 16,
    'MCP-1(phys)' : 17,
    'Insulin(ug/l)' : 18,
    'Glucose/Insulin' : 19,
    'Leptin(pg/ml)' : 20,
    'Adiponectin' : 21,
    'AorticLesions' : 22,
    'Aneurysm' : 23,
    'AorticCal-M' : 24,
    'AorticCal-L' : 25,
}



def ReadResponses():
    resps = [ [sp.strip() for sp in l.strip().split(',') ] for l in open(RESPONSE_FILE) ]
    return resps



def GetFats(resps):
    fts = []
    ids = []
    for i in range(1, len(resps)):
        mid = resps[i][RESP_RECS['MiceID']]
        tot_fat = resps[i][RESP_RECS['TotalFat']]
        if not len(tot_fat):
            ab = resps[i][RESP_RECS['AbFat']]
            ot = resps[i][RESP_RECS['OtherFat']]
            if not len(ab) or not len(ot):
                print('MISSING DATA', i, mid, 'ob', ab, 'ot', ot, '100x', resps[i][RESP_RECS['100xfat/weight']])
                continue

            tot_fat = float(ab) + float(ot)
            fts.append(tot_fat)
            ids.append(mid)
        else:
            fts.append(float(tot_fat))
            ids.append(mid)

    fats = np.array(fts)
    print('')
    print('fats stats:')
    print('len ', len(fats))
    print('min ', fats.min())
    print('max ', fats.max())
    print('mean', fats.mean())
    print('med ', np.median(fats))
    print('var ', fats.var())
    return fts, ids



def PlotFatsScatter(fats):
    plt.ylim(4, 6)
    plt.scatter(fats, [5 + 0.2 * random.uniform(0,1) - 0.1 for i in range(len(fats))])
    plt.savefig('scatter.png')
    plt.show()


def PlotFatsHistogram(fats):
    plt.ylabel('fat')
    plt.hist(fats)
    plt.savefig('fat.png')
    plt.show()



def MakeLabels(fats):
    labels = [ 0 if fats[i] < THERESHOLD else 1 for i in range(len(fats))]
    with open('LABELS.txt', 'w') as f:
        for l in labels:
            f.write('YES\n' if l else 'NO\n')



def ReadSourceNames():
    source_names = [ [sp.strip() for sp in l.strip().split(',') ] for l in open(SAMPLE_SOURCE_FILE) ]
    return source_names



def GetAccessionsFromSourceID(ids, source_names):
    acc = []
    for i in range(len(ids)):
        mid = ids[i].replace('_', '#').lower()
        try:
            midx = source_names[SOURCE_NAMES_IDX].index(mid)
            acc.append(source_names[ID_REF_IDX][midx])
        except ValueError:
            print('Could not find `' + mid + '`!')
            acc.append('')
    return acc



def DumpIDsAndAccessions(ids, accessions):
    with open('IDs2Accessions.txt', 'w') as f:
        for i in range(len(ids)):
            f.write('{},{}\n'.format(ids[i], accessions[i]))



def main():
    resps = ReadResponses()
    fats, ids = GetFats(resps)
    MakeLabels(fats)
    source_names = ReadSourceNames()
    accessions = GetAccessionsFromSourceID(ids, source_names)
    DumpIDsAndAccessions(ids, accessions)


if __name__ == '__main__':
    main()


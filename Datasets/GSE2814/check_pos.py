
def CheckMarkers():
    pos = [ l.split(',')[1] + l.split(',')[2] for l in open('marker.csv', 'r') ]
    pos = sorted(pos)
    uniq = list(set(pos))
    print('recs', len(pos), 'uniques', len(uniq))
    for i in range(len(uniq)):
        cnt = pos.count(uniq[i])
        if cnt != 1:
            print(i, uniq[i], cnt)



def CheckSNPs():
    with open('SNP-indicators.txt', 'r') as f:
        snps = f.readline().strip().split()
    uniq = list(set(snps))
    print('len', len(snps), 'len/2', len(snps) / 2, 'unique', len(uniq))



def main():
    CheckMarkers()
    CheckSNPs()


if __name__ == '__main__':
    main()


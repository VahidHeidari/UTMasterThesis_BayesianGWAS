import os
import time

import pandas as pd



SNP_START_FIELD = '1137'
SNP_LIST = [ 'A', 'C', 'G', 'T' ]

COL_MISSING_PERCENT = 0.90
ROW_MISSINGS_PERCENT = 0.05



def ExcelToCSV(path):
    csv_path = os.path.splitext(path)[0] + '.csv'
    if os.path.isfile(csv_path):
        print('The `' + path + '` is already converted!')
        return

    print(time.ctime() + ' Reading `' + path + '` file . . .')
    read_file = pd.read_excel(path)
    print(time.ctime() + ' Convert to `' + csv_path + '` . . .')
    read_file.to_csv(csv_path, index=None, header=True)
    print(time.ctime() + ' End of conversion!')



def ReadCSV(path):
    csv = [l.strip().split(',') for l in open(path, 'r')]
    return csv



def CountRecSNPs(rec, hdr):
    snp_start_idx = hdr.index(SNP_START_FIELD)
    snp_end_idx = len(hdr) - 4

    cnts = [0, 0, 0, 0]
    for i in range(snp_start_idx, snp_end_idx):
        if len(rec[i]):
            cnts[ ['A', 'C', 'G', 'T'].index(rec[i]) ] += 1
    return cnts



def CountColSNPs(csv):
    start_idx = csv[1].index(SNP_START_FIELD)
    end_idx = len(csv[1]) - 4
    cnts = [ 0 for i in range(start_idx, end_idx) ]
    for i in range(2, len(csv)):
        for j in range(start_idx, end_idx):
            if csv[i][j] in SNP_LIST:
                cnts[j - start_idx] += 1
    return cnts


def CalcMissings(cnts, csv):
    tot = len(csv) - 2
    missings = [ round(float(c) / tot, 2) for c in cnts ]
    return missings



def GetNumMissings(n, ignore_SNPs, csv):
    start_idx = csv[1].index(SNP_START_FIELD)
    num_missings = 0
    for i in range(start_idx, len(csv[1]) - 4):
        if i - start_idx in ignore_SNPs:
            continue

        if not csv[n][i] in SNP_LIST:
            num_missings += 1

    return num_missings



def SelectMissings(missings, min_percent):
    lows = []
    for i in range(len(missings)):
        if missings[i] < min_percent:
            lows.append(i)
    return lows



def main():
    ExcelToCSV('Mice and Genes.xlsx')
    ExcelToCSV('Mice and responses.xlsx')
    ExcelToCSV('mouse-marker-excel-file.xls')
    ExcelToCSV('journal.pgen.0020015.st001.XLS')
    ExcelToCSV('journal.pgen.0020015.st002.XLS')

    csv = ReadCSV('mouse-marker-excel-file.csv')
    MAX_ROW_MISSINGS = int((len(csv[1]) - 4 )* ROW_MISSINGS_PERCENT)
    print('COL_MISSING_PERCENT', COL_MISSING_PERCENT, 'MAX_ROW_MISSINGS', MAX_ROW_MISSINGS)

    # SNP missing per columns.
    cnts = CountColSNPs(csv)
    missings = CalcMissings(cnts, csv)
    low_quality_SNPs = SelectMissings(missings, COL_MISSING_PERCENT)
    ignore_indivs = []
    for i in range(2, len(csv)):
        if GetNumMissings(i, low_quality_SNPs, csv) >= MAX_ROW_MISSINGS:
            ignore_indivs.append(i)
    print('NUM IGNORED INDIVS: ', len(ignore_indivs))
    print(ignore_indivs)


if __name__ == '__main__':
    main()


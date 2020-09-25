import numpy as np
import scipy.stats as sts



def GetPValue(tbl):
    rows = len(tbl)
    cols = len(tbl[0])
    total_rows = [0 for r in xrange(rows)]
    total_cols = [0 for c in xrange(cols)]

    # Row total
    total = 0
    for r in xrange(rows):
        tot = 0
        for c in xrange(cols):
            tot += tbl[r][c]
        total += tot
        total_rows[r] = tot

    # Columns total
    for c in xrange(cols):
        tot = 0
        for r in xrange(rows):
            tot += tbl[r][c]
        total_cols[c] = tot

    sm = 0
    for r in xrange(rows):
        for c in xrange(rows):
            obs = tbl[r][c]
            exp = (total_rows[r] * total_cols[c]) / total
            sm += float((obs - exp) ** 2) / exp

    df = (rows - 1) * (cols - 1)
    p_val = 1 - sts.chi2.cdf(sm, df)
    print(sm, p_val)
    return p_val



def main():
    tbl = [
            [209, 280],
            [225, 248]
    ]

    GetPValue(tbl)

    tbl = [
            [207, 282],
            [231, 242]
    ]
    GetPValue(tbl)


if __name__ == '__main__':
    main()


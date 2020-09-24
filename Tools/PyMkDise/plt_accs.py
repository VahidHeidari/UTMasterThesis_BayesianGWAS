import os
import sys

import matplotlib.pyplot as plt



ACC_PATH = '../acc.txt'



def MakeRuns(path):
    fst_rec = None
    cnt = 1
    runs = {}
    with open(path, 'r') as f:
        for l in f:
            if not l.startswith('('):
                continue

            rec = eval(l.strip())
            if rec[0] == fst_rec:
                cnt += 1
            try:
                sm = runs[rec[0]]
                runs[rec[0]] =sm + rec[1]
            except KeyError:
                #print('NewKey:' + str(rec[0]))
                if not fst_rec:
                    fst_rec = rec[0]
                runs[rec[0]] = rec[1]

    for k in runs:
        runs[k] /= float(cnt)
    return runs



def main():
    if len(sys.argv) > 1:
        global ACC_PATH
        ACC_PATH = sys.argv[1]

    runs = MakeRuns(ACC_PATH)
    xs = sorted([k for k in runs])
    ys = [runs[k] for k in xs]
    for i in range(len(xs)):
        print('({}, {})'.format(xs[i], ys[i] * 100.0)),
    #plt.plot(xs, ys)
    #plt.savefig(os.path.splitext(ACC_PATH)[0] + '.png')
    #plt.show()



if __name__ == '__main__':
    main()


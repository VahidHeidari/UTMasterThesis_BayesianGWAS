import os
import random

import matplotlib.pyplot as plt


DATA_SET_PATH = 'data_set.txt'
FIG_PATH = 'fig.png'
NUM_DATA_SET = 100

TRUE_COEFS = [0.5, 0]
TRUE_PARAMS = [0, 0.5]
TRUE_RANGE = [0, 10]


def GetTruePos(x):
    return x * TRUE_COEFS[0] + TRUE_COEFS[1]


def MakeDataset():
    with open(DATA_SET_PATH, 'w') as f:
        for i in xrange(NUM_DATA_SET):
            x = random.uniform(TRUE_RANGE[0], TRUE_RANGE[1])
            d = random.gauss(GetTruePos(x), TRUE_PARAMS[1])
            f.write(str(x) + '\t' + str(d) + '\n')


def ReadDataset():
    xs = []
    ys = []
    with open(DATA_SET_PATH, 'r') as f:
        for l in f:
            sp = l.split('\t')
            xs.append(float(sp[0]))
            ys.append(float(sp[1].strip()))
    return xs, ys


def main():
    if not os.path.isfile(DATA_SET_PATH):
        MakeDataset()
    xs, ys = ReadDataset()

    plt.scatter(xs, ys)
    plt.plot(TRUE_RANGE, [GetTruePos(p) for p in TRUE_RANGE])
    plt.savefig(FIG_PATH)
    plt.show()


if __name__ == '__main__':
    main()

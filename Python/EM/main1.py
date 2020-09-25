import copy
import math
import os
import random

import matplotlib.pyplot as plt


DATA_SET_PATH = 'data_set.txt'
DATA_SET_LABEL_PATH = DATA_SET_PATH + '.lbl.txt'
FIG_PATH = 'fig.png'

NUM_DATA = 1000
TRUE_DISTS_WEIGHTS = [0.5, 0.5]
TRUE_DISTS_MU = [-3, 2]
TRUE_DISTS_SD = [1, 1]

MAX_ITR = 1000
ITR_LOG = 100
EPSILON = 1e-5


def MakeDataset():
    print('Making data set . . .')
    num_dists = [0 for i in xrange(len(TRUE_DISTS_WEIGHTS))]
    lbl = []
    with open(DATA_SET_PATH, 'w') as f:
        for i in xrange(NUM_DATA):
            idx = 0 if random.uniform(0, 1) < TRUE_DISTS_WEIGHTS[0] else 1
            lbl.append(idx)
            num_dists[idx] += 1
            smpl_data = random.gauss(TRUE_DISTS_MU[idx], TRUE_DISTS_SD[idx])
            f.write(str(smpl_data) + '\n')

    print('Dump labels . . .')
    with open(DATA_SET_LABEL_PATH, 'w') as f:
        for i in xrange(NUM_DATA):
            f.write(str(lbl[i]) + '\n')
    print('num_dists -> [' + ', '.join([str(n) for n in num_dists]) + ']')


def ReadDataset():
    data_set = []
    with open(DATA_SET_PATH, 'r') as f:
        for l in f:
            data_set.append(float(l))
    return data_set


def GetNormPDF(x, mu, sd):
    res = 1.0 / math.sqrt(2 * math.pi * sd * sd)
    pw = 0.5 * 1.0 / (sd * sd) * math.pow(x - mu, 2)
    res *= math.exp(-pw)
    return res


def IsConverged(old_theta, theta):
    for i in xrange(len(old_theta)):
        if abs(old_theta[i] - theta[i]) > EPSILON:
            return False
    return True


def main():
    # Read data set.
    if not os.path.isfile(DATA_SET_PATH):
        MakeDataset()
    data_set = ReadDataset()
    print('min', min(data_set), 'max', max(data_set))

    # Random initialize parameters.
    theta = [random.uniform(-10, 10) for i in xrange(len(TRUE_DISTS_MU))]
    assign = [[0 for i in xrange(len(theta))] for j in xrange(len(data_set))]
    print('initial theta -> [' + ', '.join([str(p) for p in theta]) + ']')

    # Run EM algorithm.
    for itr in xrange(MAX_ITR):
        old_theta = copy.copy(theta)
        if (itr + 1) % ITR_LOG == 0:
            print('iteration #{}   -> ['.format(itr + 1) + ', '.join([str(p) for p in theta]) + ']')

        # Expectation step, Fixed Params and calculate Theta
        sm = [0 for j in xrange(len(theta))]
        for i in xrange(len(data_set)):
            for j in xrange(len(theta)):
                prb = GetNormPDF(data_set[i], theta[j], 1)
                assign[i][j] = prb      # Assignment probability
                sm[j] += prb

        # Maximization step
        for j in xrange(len(theta)):
            theta[j] = 0
            for i in xrange(len(data_set)):
                theta[j] += assign[i][j] * data_set[i]
            theta[j] /= sm[j]

        if IsConverged(old_theta, theta):
            print('Converged at #{} iteration'.format(itr))
            break

    # Print inferred parameters.
    print('theta -> [' + ', '.join([str(round(p, 2)) for p in theta]) + ']')

    # Plot histogram.
    plt.hist(data_set, bins=20)
    plt.savefig(FIG_PATH)
    plt.show()


if __name__ == '__main__':
    main()

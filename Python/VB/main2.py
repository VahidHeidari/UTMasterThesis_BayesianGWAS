import os
import math
import random

import matplotlib.pyplot as plt


DATA_SET_PATH = 'nrm_data_set.txt'
IS_SHOW_HISTOGRAM = True
NUM_DATA_SET = 1000
TRUE_MU = 10.0
TRUE_SD = 3
NUM_ITRS = 100
EPSILON = 0.00001


def MakeDataset():
    with open(DATA_SET_PATH, 'w') as f:
        for i in xrange(NUM_DATA_SET):
            f.write(str(random.gauss(TRUE_MU, TRUE_SD)) + '\n')


def main():
    if not os.path.isfile(DATA_SET_PATH):
        MakeDataset()
    data_set = [float(l) for l in open(DATA_SET_PATH, 'r')]

    sum_xi = 0
    sum_xi_2 = 0
    for xi in data_set:
        sum_xi += xi
        sum_xi_2 += xi * xi

    N = len(data_set)
    mu_0 = 0.0
    k_0 = 1.0
    a_0 = 1.0
    b_0 = 2.0

    mu_n = ((k_0 * mu_0) + sum_xi) / (k_0 + N)
    a_n = a_0 + ((N + 1.0) / 2.0)

    #tu_n = random.gammavariate(a_0, b_0)
    tu_n = 1.0 / random.uniform(0.1, 10)
    b_n = 1
    for itr in xrange(NUM_ITRS):
        # Save previous parameters for testing convergence.
        prev_b_n = b_n
        prev_tu_n = tu_n
        print('itr #{}    b_n {:.2f}    tu_n {:}'.format(itr, b_n, tu_n))

        exp_mu = mu_n
        exp_mu_2 = (1.0 / tu_n) + math.pow(mu_n, 2)
        mu_tu = a_n / b_n
        b_n = b_0 + \
            (k_0 / 2.0) * (exp_mu_2 + math.pow(mu_0, 2) - (2 * exp_mu * mu_0)) + \
            0.5 * (sum_xi_2 + (N * exp_mu_2) - (2.0 * exp_mu * sum_xi))
        tu_n = (k_0 + N) * mu_tu

        # Test convergence.
        if abs(prev_b_n - b_n) < EPSILON and abs(prev_tu_n - tu_n) < EPSILON:
            print('Converged at itr #' + str(itr))
            break

    print('MU:')
    print('mu_n', round(mu_n, 2))
    print('tu_n', round(tu_n, 2))
    print('sd_n', round(math.sqrt(1.0 / tu_n), 2))
    print('\nSD:')
    print('a_n ', round(a_n, 2))
    print('b_n ', round(b_n, 2))
    print('sd  ', round(math.sqrt(b_n / a_n), 2))

    if IS_SHOW_HISTOGRAM:
        plt.hist(data_set, bins=10)
        plt.show()


if __name__ == '__main__':
    main()

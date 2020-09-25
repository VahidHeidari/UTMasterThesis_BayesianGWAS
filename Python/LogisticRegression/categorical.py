import os
import random

import matplotlib.pyplot as plt
import numpy as np
import sklearn.linear_model



CAT_DATA_SET_PATH = 'cat_dataset.txt'

NUM_CATS = 3
NUM_DATA = 100
STEP = 1.0 / NUM_CATS



def GetRandomCategory():
    r = random.uniform(0, 1)
    sm = 0
    for j in xrange(NUM_CATS):
        if sm + STEP >= r:
            return j
        sm += STEP
    return NUM_CATS - 1



def MakeCategoricalData():
    print('Making categorical data set.')
    with open(CAT_DATA_SET_PATH, 'w') as f:
        for i in xrange(NUM_DATA):
            k = GetRandomCategory()
            y = 1 if k == 2 else 0
            f.write('{}\t{}\n'.format(y, k))


def Predict(x, clf):
    print('PREDICT  ', x, clf.predict(x))



def main():
    if not os.path.isfile(CAT_DATA_SET_PATH):
        MakeCategoricalData()

    # Reading data set.
    data_set = [l.split() for l in open(CAT_DATA_SET_PATH, 'r')]
    y = [d[0] for d in data_set]
    x = [d[1] for d in data_set]
    bins = [[i if int(d[1]) == i else 0 for i in xrange(NUM_CATS)] for d in data_set]

    # Fit model.
    clf = sklearn.linear_model.LogisticRegression().fit(bins, y)
    print('Stats:')
    #print('clf      ', clf)
    print('score    ', clf.score(bins, y))
    #print('params   ', clf.get_params())
    print('coef     ', clf.coef_)
    print('intercept', clf.intercept_)
    print('classes  ', clf.classes_)
    print('iters    ', clf.n_iter_)
    Predict([[1, 0, 0]], clf)
    Predict([[0, 1, 0]], clf)
    Predict([[0, 0, 1]], clf)
    Predict([[1, 0, 1]], clf)
    Predict([[0, 1, 1]], clf)
    Predict([[1, 1, 1]], clf)

    #plt.scatter(x, y)
    #plt.show()


if __name__ == '__main__':
    main()


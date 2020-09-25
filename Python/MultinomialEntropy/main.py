import math

import numpy as np
import scipy.stats as st



def Fact(x):
    res = 1
    for i in range(2, x + 1):
        res *= i
    return res


def LogFact(x):
    res = 0
    for i in range(2, x + 1):
        res += np.log(i)
    return res


def Comb(n, k):
    if (k == 0 or k == n):
        return 1
    if (k == 1 or k == n - 1):
        return n

    lg = LogFact(n) - LogFact(n - k) - LogFact(k)
    return np.exp(lg)


def MultComb(n, x):
    denom = 1
    for x_i in x:
        denom *= Fact(x_i)
    res = Fact(n) / denom
    return res


def PMF(n, p, x):
    res = MultComb(n, x)
    for i in range(len(p)):
        res *= math.pow(p[i], x[i])
    return res


def Entropy(n, p):
    ent = -LogFact(n)

    sm = 0
    for i in range(len(p)):
        sm += p[i] * np.log(p[i])
    ent -= n * sm

    sm = 0
    for i in range(len(p)):
        for x_i in range(n):
            sm += Comb(n, x_i) * math.pow(p[i], x_i) * math.pow(1 - p[i], n - x_i) * LogFact(x_i)
    ent += sm
    return ent



def main():
    #n = 8
    #p = [0.3, 0.2, 0.5]
    #x = [1, 3, 4]

    #n = 8
    #p = [0.1, 0.2, 0.2, 0.3, 0.2]
    #x = [  1,   1,   2,   1,   3]

    n = 100
    p = [0.1, 0.8, 0.1]
    x = [ 10,  80,  10]

    ml = st.multinomial(n, p)
    print('ml.pmf   -> {}'.format(ml.pmf(x)))
    print('my PMF   -> {}'.format(PMF(n, p, x)))
    print('ml.entorpy -> {}'.format(ml.entropy()))
    print('my Entorpy -> {}'.format(Entropy(n, p)))


if __name__ == '__main__':
    main()


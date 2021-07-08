import math

import matplotlib.pyplot as plt


def GetSigmoid(x):
    return 1.0 / (1.0 + math.exp(-x))


def GetSigmoidDrived(x):
    s = GetSigmoid(x)
    return s * (1.0 - s)


def main():
    x = [ float(i) / 100.0 for i in range(-8 * 100, 8 * 100) ]
    y = [ GetSigmoid(x[i]) for i in range(len(x)) ]
    ys = [ GetSigmoidDrived(x[i]) for i in range(len(x)) ]
    plt.plot(x, y, x, ys)
    plt.show()


if __name__ == '__main__':
    main()


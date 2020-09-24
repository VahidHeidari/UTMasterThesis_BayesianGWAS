import math
import os
import sys

import matplotlib.pyplot as plt



ONE_MODEL = [ 220, 2392, 5505,   78, 7389, 4216]                    # One way
TWO_MODEL = [6183, 8139,  539, 3242, 5435, 9491, 2969, 8679]      # Two way
TRE_MODEL = [1724, 2620, 6005,  578, 3786, 9970]                  # Three way

K1_MODEL = [1724, 2620, 6005]
K2_MODEL = [ 578, 3786, 9970]



def PlotDiseLoci(xs, ys, dise_loci, color):
    plt.scatter([x * 0.1 for x in dise_loci], [ys[i] for i in dise_loci], color=color)
    for i in dise_loci:
        plt.annotate(str(i), (xs[i], ys[i]), color=color)


if __name__ == '__main__':
    file_path = 'p_vals.txt' if len(sys.argv) < 2 else sys.argv[1]
    print('Try opening `' + file_path + '\' file . . .')
    xs = []
    ys = []
    with open(file_path, 'r') as f:
        for l in f:
            l = l.strip()
            xs.append(float(l.split('\t')[0]))
            ys.append(float(l.split('\t')[1]))

    plt.ylim(-0.5, 11.5)
    plt.xlabel('Locus')
    plt.ylabel('-log(p-value)')
    plt.scatter(xs, ys)

    #for i in range(len(ys)):
    #    if ys[i] > -math.log(0.05):
    #        x = int(xs[i] * 10)
    #        plt.annotate(str(int(x)), (xs[i], ys[i]))

    PlotDiseLoci(xs, ys, ONE_MODEL, 'black')
    PlotDiseLoci(xs, ys, TWO_MODEL, 'orange')
    PlotDiseLoci(xs, ys, TRE_MODEL, 'red')

    #PlotDiseLoci(xs, ys, K1_MODEL, 'red')
    #PlotDiseLoci(xs, ys, K2_MODEL, 'red')

    plt.savefig(os.path.splitext(file_path)[0] + '.png')
    plt.show()


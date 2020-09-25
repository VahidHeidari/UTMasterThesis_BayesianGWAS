import math
import random
import sys

import matplotlib.pyplot as plt



JITTER = 3



def main():
    st = int(sys.argv[1])
    md = int(sys.argv[2])
    ed = int(sys.argv[3])
    xs = []
    ys = []
    with open('data.txt', 'w') as f:
        for i in range(101):
            dt = random.uniform(st, ed)
            y = 0 if dt < md else 1
            dt += random.uniform(-JITTER, JITTER)
            xs.append(dt)
            ys.append(y)
            f.write(str(round(dt, 2)) + '\t' + str(y) + '\n')
    plt.scatter(xs, ys)
    plt.show()


if __name__ == '__main__':
    main()


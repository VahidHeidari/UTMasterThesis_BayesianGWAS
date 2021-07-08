import random

import matplotlib.pyplot as plt
import numpy as np



MAX_ITR = 100000
LOG_ITR = max(1, MAX_ITR // 10)
EPSILON = 0.000001
LEARNING_RATE = 0.001

params = []



def GetSigmoid(params, xs):
    sm = params[len(params) - 1]            # Bias is the last element in params array.
    for i in range(len(xs)):
        sm += params[i] * xs[i]
    res = 1.0 / (1.0 + np.exp(-sm))
    return res


def GetSigmoidDrived(params, xs):
    s = GetSigmoid(params, xs)
    res = s * (1.0 - s)
    return res


def GetLoss(params, data):
    j_sm = 0
    for d in data:
        xs = d[0 : len(d) - 1]
        y = d[len(d) - 1]
        sm = GetSigmoid(params, xs)
        j_sm += np.log(sm) if y == 1 else np.log(1.0 - sm)
    m = float(len(data))
    j = -1.0 / m * j_sm
    return j



def ReadData(path):
    data = []
    with open(path, 'r') as f:
        for l in f:
            sp = l.split('\t')
            rec = (float(sp[0]), int(sp[1].strip()))
            data.append(rec)
    return data



def IsParamsDifferent(params, prev_params):
    for i in range(len(params)):
        if abs(params[i] - prev_params[i]) > EPSILON:
            return True

    return False



def main():
    print('Read data')
    data = ReadData('data.txt')

    print('Initilize random parameters')
    params = [ round(random.uniform(-10, 10), 2) for i in range(len(data[0])) ]
    print('    params :', params)
    print('    GetLoss:', round(GetLoss(params, data), 2))

    print('start loop . . .')
    for i in range(MAX_ITR):
        if i % LOG_ITR == 0:
            print('iteration #{}'.format(i))
            print('    params :', params)
            print('    GetLoss:', round(GetLoss(params, data), 2))


        gradian = np.zeros(len(params))
        for j in range(len(data)):
            xj = data[j][:len(data[0]) - 1]
            yj = data[j][len(data[0]) - 1]
            err = GetSigmoid(params, xj) - yj
            gradian += np.append(xj, [1]) * LEARNING_RATE * err

        new_params = params - gradian
        if not IsParamsDifferent(new_params, params):
            print('--- BREAK LOOP ---!')
            break

        params = new_params

    print('Results after #{} iterations:'.format(i))
    print('    params :', params)
    print('    GetLoss:', round(GetLoss(params, data), 2))

    tot_true = 0.0
    for i in range(len(data)):
        sig = GetSigmoid(params, data[i][:len(data[0]) - 1])
        y_pred = 1 if sig > 0.5 else 0
        y = data[i][len(data[0]) - 1]
        if y_pred == y:
            tot_true += 1.0
    acc = round(tot_true / len(data) * 100.0, 2)
    print('    Acc:', acc)

    min_range = int(min(data)[0] - 5) * 100
    max_range = int(max(data)[0] + 5) * 100
    x = [ float(i) / 100.0 for i in range(min_range, max_range) ]
    y = [ GetSigmoid(params, [x[i]]) for i in range(len(x)) ]
    plt.plot(x, y)
    x = [ d[0] for d in data ]
    y = [ d[1] for d in data ]
    plt.scatter(x, y)
    plt.show()


if __name__ == '__main__':
    main()


import copy
import os
import random

import main1


INFERRED_LABEL_PATH = 'inf_label.txt'


def CheckClusterLabels(data_set, params):
    # Read labels.
    with open(main1.DATA_SET_LABEL_PATH, 'r') as f:
        lbls = [int(l) for l in f]

    # Map labels.
    map_lbls = {0: 0, 1: 0}
    if abs(params[0] - main1.TRUE_DISTS_MU[0]) < 0.5:
        map_lbls[1] = 1
    else:
        map_lbls[0] = 1
    print('label map -> ' + str(map_lbls))

    tot_ok = 0
    with open(INFERRED_LABEL_PATH, 'w') as f:
        for i in xrange(len(data_set)):
            z = [main1.GetNormPDF(data_set[i], params[0], 1), main1.GetNormPDF(data_set[i], params[1], 1)]
            z_lbl = map_lbls[1 if z[1] >= z[0] else 0]
            f.write(str(z_lbl) + '\t' + str(z[0]) + '\t' + str(z[1]) + '\n')
            if z_lbl == lbls[i]:
                tot_ok += 1
            else:
                print('FALSE label i:' + str(i + 1) + '    z_lbl:' + str(z_lbl) + '   lbl:' + str(lbls[i]))
    print('total OK : ' + str(tot_ok) + ' of ' + str(len(data_set)))


def main():
    if not os.path.isfile(main1.DATA_SET_PATH):
        main1.MakeDataset()
    data_set = main1.ReadDataset()

    params = [random.uniform(-10, 10), random.uniform(-10, 10)]
    print('Initial params -> [' + ', '.join([str(p) for p in params]) + ']')
    z = [[0, 0] for i in xrange(len(data_set))]

    for itr in xrange(main1.MAX_ITR):
        old_params = copy.copy(params)

        # Expectation step
        for i in xrange(len(data_set)):
            z[i][0] = main1.GetNormPDF(data_set[i], params[0], 1)
            z[i][1] = main1.GetNormPDF(data_set[i], params[1], 1)

        # Maximization step
        sm_zx = [0, 0]
        sm_z = [0, 0]
        for i in xrange(len(data_set)):
            sm_zx[0] += z[i][0] * data_set[i]
            sm_zx[1] += z[i][1] * data_set[i]
            sm_z[0] += z[i][0]
            sm_z[1] += z[i][1]
        params[0] = sm_zx[0] / sm_z[0]
        params[1] = sm_zx[1] / sm_z[1]

        print('iteration #' + str(itr) + ' params -> [' + ', '.join([str(p) for p in params]) + ']')
        if main1.IsConverged(old_params, params):
            print('Converged at #{} iteration'.format(itr))
            break

    CheckClusterLabels(data_set, params)


if __name__ == '__main__':
    main()

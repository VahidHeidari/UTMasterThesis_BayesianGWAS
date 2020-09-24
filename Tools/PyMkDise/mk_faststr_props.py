import os
import sys

import mk_str_props

sys.path.append('..')
import run_faststr



DIFF_PERCENT = run_faststr.DIFF_PERCENT
INPUTS = run_faststr.INPUTS

FAST_STRUCTURE_OUTPUT_PATH = 'Outs/out_test_str.2.meanQ'
FAST_STRUCTURE_PROP_PATH = 'props.txt'



def MakeFastStructureProps(in_path, out_path):
    # Read original fast stucture output.
    with open(in_path, 'r') as f:
        qs = [l.strip().split() for l in f]

    # Calclulate cluster and write it down.
    with open(out_path, 'w') as f:
        for q in qs:
            f.write(' '.join(q))
            f.write('     K:' + str(q.index(max(q))) + '\n')



def main():
    for i in INPUTS:
        in_path = os.path.join('Outs', 'N{}_D{}_outs.2.meanQ'.format(i, DIFF_PERCENT))
        clusters = mk_str_props.ReadOrigFastStructureOutput(in_path)
        clusters_count = mk_str_props.GetClusterCounts(mk_str_props.NUM_CLUSTERS, clusters)
        accs = mk_str_props.CalcAccuracy(clusters_count)
        print(i, max(accs), accs)


if __name__ == '__main__':
    main()


import os
import subprocess
import sys
import time

sys.path.append('TestStr')
import mk_str_props



DIFF_PERCENT = 1
NUM_CLUSTERS = 3

IS_DEBUG = False
IS_KEEP_TEMPS = False

NUM_LOCI = 10000
INPUTS = [50, 100, 200, 300, 400, 500] #, 1000]

ACC_PATH = 'acc-K{}-D{}.txt'.format(NUM_CLUSTERS, DIFF_PERCENT)



def RunCommand(cmd):
    if IS_DEBUG:
        if os.name == 'nt':
            cmd = ['cmd', '/C', 'echo'] + cmd
        else:
            cmd = ['echo'] + cmd

    cmd = [str(itm) for itm in cmd]
    print(' START : ' + time.ctime(time.time()))
    print(' '.join(cmd))
    ret = subprocess.call(cmd)
    print('ret : ' + str(ret))
    print(' END   : ' + time.ctime(time.time()))
    return ret



def DeleteTemps(cpp_geno, i):
    if IS_KEEP_TEMPS:
        return

    if os.path.isfile(cpp_geno):
        os.remove(cpp_geno)

    cpp_freq = 'cpp_freqs_K{}_L{}_D{}_N{}.txt'.format(NUM_CLUSTERS, NUM_LOCI, DIFF_PERCENT, i)
    if os.path.isfile(cpp_freq):
        os.remove(cpp_freq)

    fst_genos = cpp_geno + '-faststr.txt'
    if os.path.isfile(fst_genos):
        os.remove(fst_genos)



def AppendAcc(acc_path, props_path, num_indivs):
    with open(acc_path, 'a') as f:
        clusters = mk_str_props.ReadOrigFastStructureOutput(props_path)
        clusters_count = mk_str_props.GetClusterCounts(NUM_CLUSTERS, clusters)
        accs = mk_str_props.CalcAccuracy(clusters_count)
        f.write(str((num_indivs, max(accs), accs)) + '\n')



def RunFastStructure():
    for i in INPUTS:
        print('Running ' + str(i) + ' . . .\n')

        mk_str_cmds = ['./mk_str.out', NUM_CLUSTERS, NUM_LOCI, i, DIFF_PERCENT]
        RunCommand(mk_str_cmds)

        cpp_geno = 'cpp_genos_K{}_L{}_D{}_N{}.str'.format(NUM_CLUSTERS, NUM_LOCI, DIFF_PERCENT, i)
        out_file = os.path.join('TestStr', 'Outs', 'N{}_K{}_D{}_outs'.format(i, NUM_CLUSTERS, DIFF_PERCENT))
        cmds = ['python2', 'structure.py', '-K', NUM_CLUSTERS, '--input=' + os.path.splitext(cpp_geno)[0],
                '--format=str', '--output=' + out_file, '--tol=1e-2']#, '--full']
        RunCommand(cmds)

        DeleteTemps(cpp_geno, i)
        print('----------------\n\n')
        AppendAcc(ACC_PATH, '{}.{}.meanQ'.format(out_file, NUM_CLUSTERS), i)

    with open(ACC_PATH, 'a') as f:
        f.write('\n')



def main():
    for i in range(5):
        print('------------------------------------')
        print('           ITR: ' + str(i + 1))
        print('------------------------------------')
        print('')
        RunFastStructure()


if __name__ == '__main__':
    main()


import os
import sys



BASE_PATH = 'D:\\' if os.name == 'nt' else '/media/vahidlinux20/'
sys.path.append(os.path.join(BASE_PATH, 'Datasets', 'Programs', 'Mapper'))
sys.path.append(os.path.join('..', 'fastStructure', 'TestStr'))
import mk_str_props
import run_cygwin



DIFF_PERCENT = 1
NUM_CLUSTERS = 3

IS_DEBUG = False
IS_KEEP_TEMPS = False

NUM_LOCI = 1000
INPUTS = [50, 100, 200, 300]#, 400, 500]

ACC_PATH = 'acc-K{}-D{}.txt'.format(NUM_CLUSTERS, DIFF_PERCENT)
STRUCT_EXE = os.path.join('structure_kernel_src', 'structure.exe')
MK_STR_EXE = os.path.join('..', 'fastStructure', 'mk_str' + ('.out' if os.name == 'posix' else '.exe'))

MAIN_PARAMS = '''#define MAXPOPS  {}
#define INFILE   {}
#define OUTFILE  {}
#define NUMINDS  {}
#define NUMLOCI  {}
'''



def RunCommand(cmd):
    if IS_DEBUG:
        if os.name == 'nt':
            cmd = ['cmd', '/C', 'echo'] + cmd
        else:
            cmd = ['echo'] + cmd

    cmd = [str(itm) for itm in cmd]
    return run_cygwin.RunCommand(cmd)



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
        props = mk_str_props.ReadOrigStructureOutput(props_path)
        clusters = [rec.index(max(rec)) for rec in props]
        clusters_count = mk_str_props.GetClusterCounts(NUM_CLUSTERS, clusters)
        accs = mk_str_props.CalcAccuracy(clusters_count)
        f.write(str((num_indivs, max(accs), accs)) + '\n')



def RunStructure():
    for i in INPUTS:
        cpp_geno = 'cpp_genos_K{}_L{}_D{}_N{}.str'.format(NUM_CLUSTERS, NUM_LOCI, DIFF_PERCENT, i)
        out_file = os.path.join('Outs', 'N{}_K{}_D{}_outs'.format(i, NUM_CLUSTERS, DIFF_PERCENT))

        mk_str_cmds = [MK_STR_EXE, NUM_CLUSTERS, NUM_LOCI, i, DIFF_PERCENT]
        RunCommand(mk_str_cmds)

        with open('mainparams', 'w') as f:
            f.write(MAIN_PARAMS.format(NUM_CLUSTERS, cpp_geno, out_file, i * NUM_CLUSTERS, NUM_LOCI))
        RunCommand([STRUCT_EXE])

        DeleteTemps(cpp_geno, i)
        AppendAcc(ACC_PATH, out_file + '_f', i)

    with open(ACC_PATH, 'a') as f:
        f.write('\n')



def main():
    run_cygwin.AddCygwinPath()
    for i in range(1):
        print('------------------------------------')
        print('           ITR: ' + str(i + 1))
        print('------------------------------------')
        print('')
        RunStructure()
    run_cygwin.RemoveCygwinPath()


if __name__ == '__main__':
    main()


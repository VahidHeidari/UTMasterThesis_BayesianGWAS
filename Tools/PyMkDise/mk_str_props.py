import os
import itertools



LINUX_BASE = '/media/vahidlinux20/2CB4AE50B4AE1BF8/'

# Original structure
#STRUCT_OUT_PATH_WIN = 'F:\\C++\\Structure\\structure_kernel_src\\genos_f'
STRUCT_OUT_PATH_WIN = 'F:\\C++\\Structure\\structure_kernel_src\\genos_no_ext_cols.out_f'
STRUCT_OUT_PATH_LINUX = LINUX_BASE + 'C++/Structure/GCC/structure_kernel_src/genos_f'
STRUCT_OUT_PATH = STRUCT_OUT_PATH_LINUX if os.name == 'posix' else STRUCT_OUT_PATH_WIN
PROPS_OUT_PATH = 'str_props.txt'

STR_OFF = 5
STR_OFF = 4

# Original fasta structure
#FAST_STRUCT_OUT_PATH = 'Outs/out_test_str.2.meanQ'
#FAST_STRUCT_OUT_PATH = 'Outs/cpp_genos.2.meanQ'
FAST_STRUCT_OUT_PATH = 'N50/outs.2.meanQ'

# My structures
MY_FAST_STRUCT_OUT_PATH_WIN = 'F:\\C++\\MySTRUCTURE\\props.txt'
MY_FAST_STRUCT_OUT_PATH_LINUX = LINUX_BASE + 'C++/MySTRUCTURE/props.txt'
MY_FAST_STRUCT_OUT_PATH = MY_FAST_STRUCT_OUT_PATH_WIN if os.name == 'nt' else MY_FAST_STRUCT_OUT_PATH_LINUX
PROPS_MY_FAST_STRUCT_OUT_PATH = 'my_fast_str_props.txt'
MY_NO_ADMIX_OUT_PATH = 'F:\\C++\\MySTRUCTURE\\no_admix_props.txt'

NUM_CLUSTERS = 2



def ReadOrigStructureOutput(path):
    props = []
    with open(path, 'r') as f:
        l = f.readline().strip()
        while l != 'Inferred ancestry of individuals:':
            l = f.readline().strip()

        l = f.readline()
        l = f.readline().strip()
        while len(l) != 0:
            props.append([float(q) for q in l.split()[STR_OFF:]])
            l = f.readline().strip()
    return props



def ReadOrigFastStructureOutput(path):
    clusters = []
    with open(path, 'r') as f:
        for l in f:
            qs = [float(q) for q in l.strip().split()]
            clusters.append(qs.index(max(qs)))
    return clusters



def ReadMyFastStructureOutput(path):
    props = []
    with open(path, 'r') as f:
        for l in f:
            l = l.strip()
            if len(l) == 0:
                continue

            qs = []
            for sp in l.split()[1:]:
                if sp == 'Q:':
                    break
                qs.append(float(sp))
            props.append(qs)
    return props



def WriteProps(path, props):
    clusters = []
    with open(path, 'w') as f:
        for rec in props:
            mx_k = rec.index(max(rec))
            clusters.append(mx_k)
            f.write(str(rec) + '    K:' + str(mx_k) + '\n')
    return clusters



def ReadMyNoAdmix(path):
    clusters = []
    with open(path, 'r') as f:
        l = f.readline()
        for l in f:
            l = l.strip()
            if len(l) == 0:
                break
            clusters.append(int(l.split()[1]))
    return clusters



def GetClusterCounts(num_clusters, clusters):
    clstr_count = [[0 for k in range(num_clusters)] for i in range(num_clusters)]
    num_indivs_in_cluster = len(clusters) // num_clusters 
    for k in range(num_clusters):
        for i in range(num_indivs_in_cluster):
            indiv_idx = k * num_indivs_in_cluster + i
            clstr_idx = clusters[indiv_idx]
            clstr_count[k][clstr_idx] += 1
    return clstr_count



def GetAcc(prm, cnts):
    tot = 0
    for rec in cnts:
        for c in rec:
            tot += c

    num_clstr = len(cnts[0])
    sm = 0.0
    for k in range(num_clstr):
        idx = prm[k]
        sm += cnts[k][idx]
    acc = sm / tot
    return acc



def CalcAccuracy(clusters_count):
    num_clstr = len(clusters_count[0])
    items = [j for j in range(num_clstr)]
    accs = [GetAcc(prm, clusters_count) for prm in itertools.permutations(items)]
    return accs



def main():
    # Original structure output
    #props = ReadOrigStructureOutput(STRUCT_OUT_PATH)
    #clusters = WriteProps(PROPS_OUT_PATH, props)

    # Original fast structure output
    clusters = ReadOrigFastStructureOutput(FAST_STRUCT_OUT_PATH)

    # My fast structure output
    #props = ReadMyFastStructureOutput(MY_FAST_STRUCT_OUT_PATH)
    #clusters = WriteProps(PROPS_MY_FAST_STRUCT_OUT_PATH, props)

    # My no admixture output
    #clusters = ReadMyNoAdmix(MY_NO_ADMIX_OUT_PATH)

    clusters_count = GetClusterCounts(NUM_CLUSTERS, clusters)
    accs = CalcAccuracy(clusters_count)
    print(max(accs), accs)


if __name__ == '__main__':
    main()


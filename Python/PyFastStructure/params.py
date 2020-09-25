import random

import inputs
import numpy as np


P = []
Z = []
Q = []


# ######### Z ##########
def MakeNoAdmixAssignments(data_set):
    global Z
    N = inputs.GetNumIndivs(data_set)
    Z = [int(random.uniform(0, inputs.NUM_CLUSTERS)) for _ in xrange(N)]
    return Z


def GetNoAdmixAssignment(assignment, individual_idx):
    return assignment[individual_idx]


# ---------- Z, Q ----------
def MakeNormalizeAdmixZ(assignments):
    for ind in assignments:
        for chrm in ind:
            for loc in chrm:
                sm = np.sum(loc)
                loc /= sm


def MakeAdmixAssignments(data_set):
    global Z, Q
    N = inputs.GetNumIndivs(data_set)
    L = inputs.GetNumLoci(data_set)
    C = inputs.GetNumChr(data_set)
    K = inputs.NUM_CLUSTERS
    Z = np.ones((N, C, L, K)) + (0.1 * (0.5 - np.random.rand(N, C, L, K)))
    MakeNormalizeAdmixZ(Z)
    Q = (1.0 / K) * np.ones((N, K))
    return Z, Q


def GetAdmixAssignment(assignment, individual_idx, chromosome_idx, locus_idx, cluster_idx):
    return assignment[individual_idx][chromosome_idx][locus_idx][cluster_idx]


def SetAdmixAssignment(value, assignment, individual_idx, chromosome_idx, locus_idx, cluster_idx):
    assignment[individual_idx][chromosome_idx][locus_idx][cluster_idx] = value


def GetAdmixProp(props, individual_idx, cluster_idx):
    return props[individual_idx][cluster_idx]


def SetAdmixProp(value, props, individual_idx, cluster_idx):
    props[individual_idx][cluster_idx] = value


# ---------- P ----------
def MakeFreqs(data_set):
    global P
    L = inputs.GetNumLoci(data_set)
    K = inputs.NUM_CLUSTERS
    P = 0.5 * np.ones((K, L, 2))
    return P


def GetFreq(freqs, cluster_idx, locus_idx, param_idx):
    return freqs[cluster_idx][locus_idx][param_idx]


def SetFreq(value, freqs, cluster_idx, locus_idx, params_idx):
    freqs[cluster_idx][locus_idx][params_idx] = value


if __name__ == '__main__':
    print('This is a library :(')

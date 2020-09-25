import copy

import numpy as np
import scipy.special as sp

import inputs
import params


MAX_ITR = 500
EPSILON = 1e-3


def PrintDatasets(data_set, genotypes, assignments):
    print('NumIndivs:', inputs.GetNumIndivs(data_set))
    print('NumChr:   ', inputs.GetNumChr(data_set))
    print('NumLoci:  ', inputs.GetNumLoci(data_set))
    print('NumAsgs:  ', len(assignments))
    ind = 5; loc = 8
    print('indiv:' + str(ind) + ', chr:0, locus:' + str(loc) + '  allele = ' + str(inputs.GetAllele(data_set, ind, 0, loc)))
    print('indiv:' + str(ind) + ', chr:1, locus:' + str(loc) + '  allele = ' + str(inputs.GetAllele(data_set, ind, 1, loc)))
    print('indiv:' + str(ind) + '  locus:' + str(loc) + '         genotype=' + str(params.GetGenotype(genotypes, ind, loc)))


def IsConverged(old_assignments, old_props, old_freqs, assignments, props, freqs):
    for n in xrange(len(old_assignments)):
        for c in xrange(len(old_assignments[0])):
            for l in xrange(len(old_assignments[0][0])):
                for k in xrange(len(old_assignments[0][0][0])):
                    if abs(params.GetAdmixAssignment(old_assignments, n, c, l, k) -
                           params.GetAdmixAssignment(assignments, n, c, l, k)) > EPSILON:
                        return False

    for n in xrange(len(old_props)):
        for k in xrange(len(old_props[0])):
            if abs(params.GetAdmixProp(old_props, n, k) -
                   params.GetAdmixProp(props, n, k)) > EPSILON:
                return False

    for k in xrange(len(old_freqs)):
        for l in xrange(len(old_freqs[0])):
            for i in xrange(2):
                if abs(params.GetFreq(old_freqs, k, l, i) - params.GetFreq(freqs, k, l, i)) > EPSILON:
                    return False

    return True


# Update Za_nlk, Zb_nkl.
def UpdateZ(genotypes, assignments, props, freqs):
    for n in xrange(len(assignments)):
        for l in xrange(len(assignments[0][0])):
            q_0 = np.sum(props[n])
            for k in xrange(len(assignments[0][0][0])):
                geno = inputs.GetGenotype(genotypes, n, l)
                p_u = params.GetFreq(freqs, k, l, 0)
                p_v = params.GetFreq(freqs, k, l, 1)
                q_n = params.GetAdmixProp(props, n, k)

                psi_a = (sp.digamma(p_v) if geno == 0 else 0) + \
                        (sp.digamma(p_u) if geno == 1 else 0) + \
                        (sp.digamma(p_u) if geno == 2 else 0)

                psi_b = (sp.digamma(p_v) if geno == 0 else 0) + \
                        (sp.digamma(p_v) if geno == 1 else 0) + \
                        (sp.digamma(p_u) if geno == 2 else 0)

                z_a = np.exp(psi_a - sp.digamma(p_u + p_v) + sp.digamma(q_n) - sp.digamma(q_0))
                z_b = np.exp(psi_b - sp.digamma(p_u + p_v) + sp.digamma(q_n) - sp.digamma(q_0))
                params.SetAdmixAssignment(z_a, assignments, n, 0, l, k)
                params.SetAdmixAssignment(z_b, assignments, n, 1, l, k)
    params.MakeNormalizeAdmixZ(assignments)


# Update Q.
def UpdateQ(assignments, props):
    L = len(assignments[0][0])
    K = len(props[0])
    for n in xrange(len(props)):
        for k in xrange(len(props[0])):
            sm_za_zb = 0
            for l in xrange(L):
                sm_za_zb += params.GetAdmixAssignment(assignments, n, 0, l, k) + \
                    params.GetAdmixAssignment(assignments, n, 1, l, k)
            q_nk = (1.0 / K) + sm_za_zb
            params.SetAdmixProp(q_nk, props, n, k)


# Update Pu, Pv.
def UpdateP(genotypes, assignments, freqs):
    # Priors
    beta = 1
    gamma = 1

    N = len(assignments)
    for k in xrange(len(freqs)):
        for l in xrange(len(freqs[0])):
            sm_za = 0
            sm_zb = 0
            for n in xrange(N):
                geno = inputs.GetGenotype(genotypes, n, l)
                z_a = params.GetAdmixAssignment(assignments, n, 0, l, k)
                z_b = params.GetAdmixAssignment(assignments, n, 1, l, k)
                sm_za += (z_a if geno == 1 else 0) + (z_a + z_b if geno == 2 else 0)
                sm_zb += (z_b if geno == 1 else 0) + (z_a + z_b if geno == 0 else 0)
            p_u = beta + sm_za
            p_v = gamma + sm_zb
            params.SetFreq(p_u, freqs, k, l, 0)
            params.SetFreq(p_v, freqs, k, l, 1)


def WriteToFile(assignments, props, freqs):
    # Marginal parameters:
    with open('Z.txt', 'w') as f:
        for n in xrange(len(assignments)):
            f.write(str(assignments[n]) + '\n\n')
    with open('Q.txt', 'w') as f:
        f.write(str(props))
    with open('P.txt', 'w') as f:
        f.write(str(freqs))

    # Assignment proportions.
    N = len(props)
    with open('Props.txt', 'w') as f:
        for n in xrange(N):
            qn = props[n]
            q0 = np.sum(qn)
            exps = ['{:0.2f}'.format(round(p / q0, 2)) for p in qn]
            mxf = max(exps); mxk = exps.index(mxf)
            f.write(', '.join(exps) + '       K:' + str(mxk) + '\n')


def main():
    #data_set = inputs.ReadGenotypes()
    data_set = inputs.ReadBinaryGenotypes()
    genotypes = inputs.MakeGenotypes(data_set)
    assignments, props = params.MakeAdmixAssignments(data_set)
    freqs = params.MakeFreqs(data_set)

    for itr in xrange(MAX_ITR):
        print('itr #' + str(itr + 1))
        old_assignments = copy.copy(assignments)
        old_props = copy.copy(props)
        old_freqs = copy.copy(freqs)

        UpdateP(genotypes, assignments, freqs)
        UpdateZ(genotypes, assignments, props, freqs)
        UpdateQ(assignments, props)

        if IsConverged(old_assignments, old_props, old_freqs, assignments, props, freqs):
            print('Converged!')
            break

    WriteToFile(assignments, props, freqs)


if __name__ == '__main__':
    main()

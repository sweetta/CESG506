import numpy as np


def get_strain(l, L):
    stretchSquared = np.dot(l,l)/np.dot(L,L)
    strain = 0.5*np.log(stretchSquared)
    return strain

def get_n(l):
    length = np.linalg.norm(l)
    n = l/length
    return n

def get_ke(EA, l, strain):
    dim = l.size
    f_int = EA*strain
    length = np.linalg.norm(l)
    n = l/length
    ke = (EA/length)*np.outer(n, n) + \
         (f_int/length)*(np.identity(dim) - np.outer(n, n))
    return ke
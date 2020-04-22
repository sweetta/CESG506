# CESG 506 - Nonlinear Analysis of Structural Systems
# Helper functions used for Main asociated with Problem 1-2
# Tatsuhiko Sweet, 4/15/2020

import numpy as np

# Returns Henkey strain from deformed and undeformed length vectors
def get_strain(l, L):
    stretchSquared = np.dot(l,l)/np.dot(L,L)
    strain = 0.5*np.log(stretchSquared)
    return strain

# Returns normal vector
def get_n(l):
    length = np.linalg.norm(l)
    n = l/length
    return n

# Returns ke of an element as derived in Problem 1-2 part 2
def get_ke(EA, l, strain):
    dim = l.size
    f_int = EA*strain
    length = np.linalg.norm(l)
    n = l/length
    ke = (EA/length)*np.outer(n, n) + \
         (f_int/length)*(np.identity(dim) - np.outer(n, n))
    return ke

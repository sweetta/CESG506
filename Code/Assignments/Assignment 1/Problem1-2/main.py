# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 1 - Problem 1-2:
# Two degree of freedom (2 DOF) problem - Load Control
# Tatsuhiko Sweet, 4/14/2020

import numpy as np
import matplotlib.pyplot as plt
from helperFunctions import *

w1 = 5.5 # m
w2 = 4.0 # m
H  = 0.5 # m
EA = 2100 # kN
max_iteration = 1000
tol = 1e-12

# Undeformed Length Vectors
L1 = np.array([w1, H])
L2 = np.array([-w2, H])

Pcr = np.array([0, -0.98171345])
gamma = [0, 0.25, 0.5, 0.75, 0.99, 0.999, 0.99999, 1]
# gamma = np.linspace(0.999999, 1.00000001, 2000)  # Used for finding Pcr

# Collect converged points as a list of tuples
# each tuple is ([Px, Py], [u, v], number of iterations, [Residual for each iteration step])
results = []
u = np.array([0, 0])
print('Tolerance = {}, s.t. tol >= |R|'.format(tol))

for g in gamma:
    P = g * Pcr
    R_list = []
    for count in range(max_iteration+1):
        l1 = L1 + u
        l2 = L2 + u
        k1 = get_ke(EA, l1, get_strain(l1, L1))
        k2 = get_ke(EA, l2, get_strain(l2, L2))
        Kf = k1 + k2
        f1 = EA*get_strain(l1, L1)*get_n(l1)
        f2 = EA*get_strain(l2, L2)*get_n(l2)
        R = P - f1 - f2
        R_list.append(np.linalg.norm(R))
        if np.linalg.norm(R) <= tol:
            print('u = {}m, (x,y) = {}m, P = {:.9f} kN = {}Pcr, '
                  'Iteration Count = {}'.format(u, L1 + u, P[1], g, count))
            results.append((P, u, count, R_list))
            break
        if count == max_iteration:
            print('u = {}m, (x,y) = {}m, P = {} kN = {}Pcr, '
                  'Stopped at Max Iteration = {}'.format(u, L1 + u, P[1], g, count))
            results.append((P, u, count, R_list, g))
        u = u + np.dot(np.linalg.inv(Kf), R)

pass
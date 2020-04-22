# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 2 - Problem 2-1:
# Displacement control for a 2 DOF problem
# Tatsuhiko Sweet, 4/21/2020

import numpy as np
import matplotlib.pyplot as plt
from helperFunctions import *

# Define geometry and material
w1 = 5.5 # m
w2 = 4.0 # m
H  = 0.5 # m
EA = 2100 # kN

# Convergence criteria
max_iteration = 10
tol = 1e-12 # kN

P0   = np.array([0, 0])
Pbar = np.array([0, -0.99])
u    = np.array([0, 0])
uk   = -np.linspace(0, 1.3, 20)
ek = np.array([0,  1])
gamma = 0

# End Inputs ##########################################################################################################

# Undeformed Length Vectors
L1 = np.array([w1, H])
L2 = np.array([-w2, H])

# Results are collect as a list of dictionaries
results = []
#   each dictionary contains:
#   'P' = P vector,
#   'u' = converged displacement vector,
#

print('Tolerance = {}, s.t. tol >= |R|'.format(tol))

step = 1

R_list = []
step_list = []
Uu = []
Uv = []
G = []

for ui in uk:
    for count in range(max_iteration+1):
        l1 = L1 + u                     # Compute Deformed Lengths
        l2 = L2 + u
        e1 = get_strain(l1, L1)         # Get Strains
        e2 = get_strain(l2, L2)
        k1 = get_ke(EA, l1, e1)         # Get ke associated
        k2 = get_ke(EA, l2, e2)
        Kf = k1 + k2                    # Get stiffness at free node
        f1 = EA*e1*get_n(l1)            # Find element internal force vectors
        f2 = EA*e2*get_n(l2)
        P  = P0 + gamma*Pbar
        R = P - f1 - f2                 # Compute residual vector
        g = np.dot(u, ek)-ui
        Rg = np.append([R], [g])
        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-ek, 0)))

        if np.linalg.norm(Rg) <= tol:        # Check for convergence using tol >= |R|
            print('u = {}m, P = {:.9f} kN = {}Pbar, '
                  'Iteration Count = {}'.format(u, P[1], gamma, count))
            results.append({'P': P, 'u': u, 'Gamma': gamma})
            Uu.append(u[0])
            Uv.append(u[1])
            G.append(gamma)
            break

        if count == max_iteration:          # Exit loop if max iteration is reached before convergence
            print('u = {}m, P = {} kN = {}Pbar, '
                  'Stopped at Max Iteration = {}'.format(u, P[1], gamma, count))
            results.append({'P': P, 'u': u, 'Gamma': gamma})

        dU = np.dot(np.linalg.inv(Kg), Rg)
        u = u + dU[[0, 1]]
        gamma = gamma + dU[2]


# Plot of gamma vs displacement, part 3 of assignment problem
plt.plot(Uu, G, label="$U_u$", marker='.')
plt.plot(Uv, G, label="$U_v$", marker='.')

plt.xlabel('displacement [m]')
plt.ylabel('$\gamma$')

plt.legend(loc='upper center', ncol=1, framealpha=1)
# plt.xlim([0, umax])
# plt.ylim([-1.5, 2])
plt.axhline(y=0, color='black')
plt.grid(True)
# plt.savefig("CESG506/Code/Assignments/Assignment 1/Problem1-1/{}".format(fileName))
plt.show()

pass
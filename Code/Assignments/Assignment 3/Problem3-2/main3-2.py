# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 3 - Problem 3-2: 2 DOF Problem - Arc-length Method
# Tatsuhiko Sweet, 4/28/2020

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
sys.path.insert(0, 'CESG506/Code/Assignments/Assignment 3')
from trussClass import *
from trussAssemble import *


Pref = 0.99    # Reference Load [kN]

# Convergence criteria
a = 0.0
max_iteration = 50
tol = 5e-11     # np.linalg.norm([R],[g]) <= tol

# Create Nodes
X = [np.array([0.0, 0.0]), np.array([5.5, 0.5]), np.array([9.5, 0.0])]

fixed_dof = [0, 1, 4, 5]
free_dof = [2, 3]

dim = len(X[0])
ndof = dim*len(X)

# Create Trusses
T = [trussClass([0, 1], X[0], X[1], 2100), trussClass([2, 1], X[2], X[1], 2100)]

# P vector constructed from reference load: P = P0 + gamma * Pbar
P0   = np.zeros(ndof)
Pbar = np.zeros(ndof)
Pbar[dim*1+1] = -Pref


P0 = P0[free_dof]
Pbar = Pbar[free_dof]

# Start at u = zero, gamma = 0.2
u_n = np.zeros(ndof)
gamma_n = 0.0
gamma_guess = 0.2

# Initialization Lists for plotting
zeroU = np.zeros(ndof)
u_guess = np.zeros(ndof)
Results = [{'u': u_n, 'gamma': gamma_n, 'R_list': 0, 'Iterations': 0}]
u_list = np.array(u_n)
Uu = [0]
Uv = [0]
G = [0]

# Initial Step to set delta_s^2
# Set Displacements and Construct Global K Matrix
[K_global, Fint] = trussAssemble(T, u_n, dim, ndof)
Kf = K_global[free_dof, :][:, free_dof]

P = P0 + gamma_guess * Pbar
R = P - Fint[free_dof]

u_guess[free_dof] = u_n[free_dof] + np.dot(np.linalg.inv(Kf), R)

s2 = np.dot(u_guess[free_dof], u_guess[free_dof]) + a*gamma_guess*gamma_guess
print('(delta_s)^2 = {}, delta_s = {}'.format(s2, np.sqrt(s2)))

flag = 0
for i in range(60):
    for count in range(max_iteration + 1):
        [K_global, Fint] = trussAssemble(T, u_guess, dim, ndof)

        # Construct Kg Matrix
        Kf = K_global[free_dof, :][:, free_dof]
        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-2*(u_guess-u_n)[free_dof], -2*a*(gamma_guess-gamma_n))))

        # Construct Rg Vector
        P = P0 + gamma_guess * Pbar
        R = P - Fint[free_dof]
        g = np.dot((u_guess-u_n)[free_dof], (u_guess-u_n)[free_dof]) + a*(gamma_guess-gamma_n)*(gamma_guess-gamma_n) - s2
        Rg = np.append([R], [g])

        # Solve for Change in u and gamma guess
        dU = np.dot(np.linalg.inv(Kg), Rg)

        # Update Guess
        u_guess[free_dof] = u_guess[free_dof] + dU[0:len(free_dof)]
        gamma_guess = gamma_guess + dU[-1]

        # Check Convergence
        if np.linalg.norm(Rg) <= tol:
            print('u = {}m, P = {:.6f} kN = {:.6f} Pbar, '
                  'Iteration Count = {}, u_{}'.format(u_n, gamma_n*Pref, gamma_n, count, i))
            break

        # Save and Print Results if it does not converge by max iteration
        if count == max_iteration:  # Exit loop if max iteration is reached before convergence
            print('u = {}m, P = {:.6f} kN = {:.6f} Pbar, '
                  'Stopped at Max Iteration = {}'.format(u_guess, gamma_n*Pref, gamma_guess, count))
            flag = 1

    if flag == 1:
        break

    # Store Converged Results
    u_list = np.vstack((u_list, u_guess))
    Uu.append(u_guess[2])
    Uv.append(u_guess[3])
    G.append(gamma_guess)

    # Initiate Next Guess
    u_last = u_n
    u_n = u_guess
    u_guess = 2*u_n - u_last
    gamma_last = gamma_n
    gamma_n = gamma_guess
    gamma_guess = 2*gamma_n - gamma_last

# Transpose for easier plotting
u_list = u_list.transpose()

# ----------------------------------------------------------------------------------------------------------------------
# Plot of gamma vs displacement, part 4 of assignment problem
# ----------------------------------------------------------------------------------------------------------------------

# Read results from HW2
HW2G = []
HW2Uu = []
HW2Uv = []
f = open("CESG506/Code/Assignments/Assignment 3/Problem3-2/FromHW2/HW2outputs", "r")
for line in f:
    row = line.strip('\n').split(',')
    HW2G.append(float(row[0]))
    HW2Uu.append(float(row[1]))
    HW2Uv.append(float(row[2]))
f.close()

plt.figure(1, figsize=(16, 8))
plt.plot(HW2Uu, HW2G, label="$U_u$ - Displacement Control", marker='o', linestyle='--', linewidth='1')
plt.plot(HW2Uv, HW2G, label="$U_v$ - Displacement Control", marker='o', linestyle='--', linewidth='1')
plt.plot(Uu, G, label="$U_u$ - Arc-length Method", marker='x', linestyle='')
plt.plot(Uv, G, label="$U_v$ - Arc-length Method", marker='x', linestyle='')

plt.xlabel('displacement [m]')
plt.ylabel('$\gamma$')

plt.legend(loc='best', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)
plt.savefig("CESG506/Code/Assignments/Assignment 3/Problem3-2/Prob3-2_GammaVsDisp")
plt.show()

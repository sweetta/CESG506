# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 2 - Problem 2-1: Using Truss Class
# Displacement control for a 2 DOF problem
# Tatsuhiko Sweet, 4/21/2020

import numpy as np
import matplotlib.pyplot as plt
from trussClass import *


X1 = np.array([0.0, 0.0])
X2 = np.array([5.5, 0.5])
X3 = np.array([9.5, 0.0])

t1 = trussClass([1, 2], X1, X2, 2100)
t2 = trussClass([3, 2], X3, X2, 2100)

elements = [t1, t2]

# Convergence criteria
max_iteration = 10
tol = 1e-12 # kN

P0   = np.array([0, 0])
Pbar = np.array([0, -0.99])
u    = np.array([0, 0])
uk   = -np.linspace(0, 1.3, 20)
ek = np.array([0,  1])
gamma = 0
zeroU = np.array([0, 0])
Uu = []
Uv = []
G = []

for ui in uk:
    for count in range(max_iteration+1):
        for truss in elements:
            truss.setDisp(zeroU, u)
        Kf = t1.get_ke() + t2.get_ke()                    # Get stiffness at free node
        P  = P0 + gamma*Pbar
        R = P - t1.get_F() - t2.get_F()                 # Compute residual vector
        g = np.dot(u, ek)-ui
        Rg = np.append([R], [g])
        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-ek, 0)))

        if np.linalg.norm(Rg) <= tol:        # Check for convergence using tol >= |R|
            print('u = {}m, P = {:.9f} kN = {}Pbar, '
                  'Iteration Count = {}'.format(u, P[1], gamma, count))
            Uu.append(u[0])
            Uv.append(u[1])
            G.append(gamma)
            break

        if count == max_iteration:          # Exit loop if max iteration is reached before convergence
            print('u = {}m, P = {} kN = {}Pbar, '
                  'Stopped at Max Iteration = {}'.format(u, P[1], gamma, count))

        dU = np.dot(np.linalg.inv(Kg), Rg)
        u = u + dU[[0, 1]]
        gamma = gamma + dU[2]
        pass

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
# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 2 - Problem 2-1:
# Going into Higher Dimensions
# Tatsuhiko Sweet, 4/21/2020

import numpy as np
import matplotlib.pyplot as plt
from trussClass import *

# ----> x
# |
# V  y          z is into page

# Nodes
X1 = np.array([0.0, 0.0, 0.0])
X2 = np.array([9.5, 0.0, 0.0])
X3 = np.array([0.0, 5.0, 0.0])
X4 = np.array([9.5, 5.0, 0.0])
X5 = np.array([5.5, 1.25, -0.5])
X6 = np.array([5.5, 3.75, -0.5])

EA = 2100

t1 = trussClass([1, 5], X1, X5, EA)
t2 = trussClass([1, 6], X1, X6, EA)
t3 = trussClass([2, 5], X2, X5, EA)
t4 = trussClass([3, 6], X3, X6, EA)
t5 = trussClass([4, 5], X4, X5, EA)
t6 = trussClass([4, 6], X4, X6, EA)
t7 = trussClass([5, 6], X5, X6, EA)

# Convergence criteria
max_iteration = 100
tol = 1e-12 # kN

P0   = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
Pbar = np.array([0.0, 0.0, 0.99, 0.0, 0.0, 0.0])
u    = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
uk   = np.linspace(0, 1.5, 200)
ek = np.array([0, 0, 1, 0, 0, 0])
gamma = 0
zeroU = np.array([0, 0, 0])
U5u = []
U5v = []
U5w = []
U6u = []
U6v = []
U6w = []
G = []

for ui in uk:
    if( (t4.get_f() > 0) and (t6.get_f() > 0)  ):
        break
    for count in range(max_iteration+1):
        u5 = u[0:3]
        u6 = u[3:6]
        t1.setDisp(zeroU, u[0:3])
        t3.setDisp(zeroU, u[0:3])
        t5.setDisp(zeroU, u[0:3])
        t2.setDisp(zeroU, u[3:6])
        t4.setDisp(zeroU, u[3:6])
        t6.setDisp(zeroU, u[3:6])
        t7.setDisp(u[0:3], u[3:6])

        K55 = t1.get_ke() + t3.get_ke() + t5.get_ke() + t7.get_ke()
        K66 = t2.get_ke() + t4.get_ke() + t6.get_ke() + t7.get_ke()
        K56 = -t7.get_ke()
        Kf = np.vstack((np.hstack((K55, K56)), np.hstack((K56, K66))))
        P  = P0 + gamma*Pbar
        R5 = P[0:3] - t1.get_F() - t3.get_F() - t5.get_F() + t7.get_F()
        R6 = -t2.get_F() - t4.get_F() - t6.get_F() - t7.get_F()
        R = np.append([R5], [R6])
        g = np.dot(u, ek)-ui
        Rg = np.append([R], [g])
        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-ek, 0)))

        if np.linalg.norm(Rg) <= tol:        # Check for convergence using tol >= |R|
            print('u = {}m, P = {:.5f} kN = {:.4f}Pbar, '
                  'Iteration Count = {}'.format(u, P[2], gamma, count))
            U5u.append(u[0])
            U5v.append(u[1])
            U5w.append(u[2])
            U6u.append(u[3])
            U6v.append(u[4])
            U6w.append(u[5])
            G.append(gamma)
            break

        if count == max_iteration:          # Exit loop if max iteration is reached before convergence
            print('u = {}m, P = {} kN = {}Pbar, '
                  'Stopped at Max Iteration = {}'.format(u, P[2], gamma, count))
        INV = np.linalg.inv(Kg)
        dU = np.dot(np.linalg.inv(Kg), Rg)
        u = u + dU[0:6]
        gamma = gamma + dU[6]
        pass

# Plot of gamma vs displacement, part 3 of assignment problem
plt.subplot(2, 1, 1)
plt.plot(U5u, G, label="$U5_u$", marker='.')
plt.plot(U5v, G, label="$U5_v$", marker='.')
plt.plot(U5w, G, label="$U5_w$", marker='.')
plt.ylabel('$\gamma$')
plt.legend(loc='best', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(U6u, G, label="$U6_u$", marker='.')
plt.plot(U6v, G, label="$U6_v$", marker='.')
plt.plot(U6w, G, label="$U6_w$", marker='.')
plt.xlabel('displacement [m]')
plt.ylabel('$\gamma$')
plt.legend(loc='best', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)

# plt.savefig("CESG506/Code/Assignments/Assignment 1/Problem1-1/{}".format(fileName))
plt.show()
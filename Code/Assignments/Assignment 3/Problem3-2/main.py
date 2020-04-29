# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 3 - Problem 3-2: 2 DOF Problem - Arc-length Method
# Tatsuhiko Sweet, 4/28/2020

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
sys.path.insert(0, 'CESG506/Code/Assignments/Assignment 3')
from trussClass import *

# Create Nodes
X1 = np.array([0.0, 0.0])
X2 = np.array([5.5, 0.5])
X3 = np.array([9.5, 0.0])

# Create Trusses
t1 = trussClass([1, 2], X1, X2, 2100)
t2 = trussClass([3, 2], X3, X2, 2100)

# Convergence criteria
a = 0.0
max_iteration = 20
tol = 1e-12 # kN

# P vector constructed from reference load: P = P0 + gamma * Pbar
P0   = np.array([0, 0])
Pbar = np.array([0, -0.99])

# Start at u = zero, gamma = 0.2
u_n = np.array([0, 0])
gamma_n = 0.0
gamma_guess = 0.2

# Initialization Lists for plotting
zeroU = np.array([0, 0])
u_list = [u_n]
Uu = [0]
Uv = [0]
G = [0]

# Initial Step to find s2

# Manually pass displacements to update Strains, n, and ke
t1.setDisp(zeroU, u_n)
t2.setDisp(zeroU, u_n)

# Manually construct Kff
Kf = t1.get_ke() + t2.get_ke()

P = P0 + gamma_guess * Pbar
R = P - t1.get_F() - t2.get_F()

u_guess = u_n + np.dot(np.linalg.inv(Kf), R)
s2 = np.dot(u_guess, u_guess) + a*gamma_guess*gamma_guess
print('(delta_s)^2 = {}, delta_s = {}'.format(s2, np.sqrt(s2)))
for i in range(100):
    for count in range(max_iteration+1):
        t1.setDisp(zeroU, u_guess)
        t2.setDisp(zeroU, u_guess)
        Kf = t1.get_ke() + t2.get_ke()

        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-2*(u_guess-u_n), -2*a*(gamma_guess - gamma_n))))

        P = P0 + gamma_guess * Pbar
        Fint = t1.get_F() + t2.get_F()
        R = P - Fint
        g = np.dot(u_guess - u_n, u_guess - u_n) + a*(gamma_guess-gamma_n)*(gamma_guess-gamma_n) - s2
        Rg = np.append([R], [g])
        dU = np.dot(np.linalg.inv(Kg), Rg)

        u_guess = u_guess + dU[[0, 1]]
        gamma_guess = gamma_guess + dU[2]
        if np.linalg.norm(Rg) <= tol:
            print('u = {}m, P = {:.9f} kN = {}Pbar, '
                              'Iteration Count = {}'.format(u_n, P[1], gamma_n, count))
            break

        # Save and Print Results if it does not converge by max iteration
        if count == max_iteration:          # Exit loop if max iteration is reached before convergence
            print('u = {}m, P = {} kN = {}Pbar, '
                  'Stopped at Max Iteration = {}'.format(u_guess, P[1], gamma_guess, count))
            Uu.append(u_guess[0])
            Uv.append(u_guess[1])
            G.append(gamma_guess)

    # Append results and create new guess before next iteration
    u_list.append(u_guess)
    Uu.append(u_guess[0])
    Uv.append(u_guess[1])
    G.append(gamma_guess)
    u_last = u_n
    u_n = u_guess
    u_guess = 2*u_n - u_last
    gamma_last = gamma_n
    gamma_n = gamma_guess
    gamma_guess = 2*gamma_n - gamma_last

    if gamma_n >= 2:
        break

# ----------------------------------------------------------------------------------------------------------------------
# Plot of gamma vs displacement, part 3 of assignment problem
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1, figsize=(16, 8))
plt.plot(Uu, G, label="$U_u$", marker='.')
plt.plot(Uv, G, label="$U_v$", marker='.')

plt.xlabel('displacement [m]')
plt.ylabel('$\gamma$')

plt.legend(loc='upper center', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)
plt.savefig("CESG506/Code/Assignments/Assignment 3/Problem3-2/Prob3-2_GammaVsDisp")
plt.show()

# # ----------------------------------------------------------------------------------------------------------------------
# # Plot of u vs v, with Fx(u,v) and Fy(u,v) contour plot, part 4 of assignment problem
# # ----------------------------------------------------------------------------------------------------------------------
# # Compute contour plot values
# grid_n = 100
# xlist = np.linspace(0.015, -0.010, grid_n)
# ylist = np.linspace(0.1, -1.3, grid_n)
# X, Y = np.meshgrid(xlist, ylist)
# Zx = np.zeros((grid_n, grid_n))
# Zy = np.zeros((grid_n, grid_n))
#
# x_idx = 0
# for xi in xlist:
#     y_idx = 0
#     for yi in ylist:
#         uv = np.array([xi, yi])
#         t1.setDisp(zeroU, uv)
#         t2.setDisp(zeroU, uv)
#         F1 = t1.get_F()
#         F2 = t2.get_F()
#         Zx[x_idx, y_idx] = F1[0] + F2[0]
#         Zy[x_idx, y_idx] = F1[1] + F2[1]
#         y_idx = y_idx + 1
#     x_idx = x_idx + 1
#
# # Plotting
# plt.figure(2, figsize=(16, 8))
#
# # Fx
# plt.subplot(1, 2, 1)
# norm_x = colors.Normalize(vmin=-20, vmax=20)
# plt.contourf(X, Y, Zx.transpose(), 100, cmap='RdBu_r', norm=norm_x)
# plt.colorbar()
# plt.plot(Uu, Uv, label='Displacement Path', marker='o', color='black')
# plt.title('u vs v with $F_x$ [kN] Overlay')
# plt.xlabel('x-displacement [m]')
# plt.ylabel('y-displacement [m]')
# plt.legend(loc='best', ncol=1, framealpha=1)
# plt.grid(True)
#
# # Fy
# plt.subplot(1, 2, 2)
# norm_y = colors.Normalize(vmin=-3, vmax=3)
# plt.contourf(X, Y, Zy.transpose(), 100, cmap='RdBu_r', norm=norm_y)
# plt.colorbar()
# plt.plot(Uu, Uv, label='Displacement Path', marker='o', color='black')
# plt.title('u vs v with $F_y$ [kN] Overlay')
# plt.xlabel('x-displacement [m]')
# plt.ylabel('y-displacement [m]')
# plt.legend(loc='best', ncol=1, framealpha=1)
# plt.grid(True)
# plt.savefig("CESG506/Code/Assignments/Assignment 2/Problem2-1/Prob2-1_ContourPlot")
#
# plt.show()
# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 3 - Problem 3-2: 2 DOF Problem - Arc-length Method
# Tatsuhiko Sweet, 4/28/2020

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
sys.path.insert(0, 'CESG506/Code/Assignments/Assignment 3')
from trussClass import *


# Structure Parameters
H = 5       # Total Height [m]
W = H/20    # Building Width [m]
dH = H/11   # Story Height [m]
EA1 = 2000  # Vertical Members [kN]
EA2 = 5000  # Horizontal Ties [kN]
EA3 = 5000  # Diagonals [kN]
Pref = 6.1685    # Reference Load [kN]

# Convergence criteria
a = 0.0
max_iteration = 50
tol = 5e-11     # np.linalg.norm([R],[g]) <= tol

# Create Nodes
X = []
for idx in range(12):
    xl = np.array([0.0, idx*dH])
    xr = np.array([W, idx*dH])
    X.append(xl)
    X.append(xr)

dim = len(X[0])
ndof = dim*len(X)

fixed_dof = [0, 1, 2, 3]
free_dof = range(4, ndof)

# Create Trusses
T = []
for idx in range(0, 21, 2):
    T.append(trussClass(  [idx, idx+2],   X[idx], X[idx+2], EA1))  # Left Vertical
    T.append(trussClass([idx+1, idx+3], X[idx+1], X[idx+3], EA1))  # Right Vertical
    T.append(trussClass(  [idx, idx+3],   X[idx], X[idx+3], EA3))  # Diagonal
    T.append(trussClass([idx+2, idx+3], X[idx+2], X[idx+3], EA2))  # Horizontal tie

del xl, xr, idx
del EA1, EA2, EA3, H

# P vector constructed from reference load: P = P0 + gamma * Pbar
P0   = np.zeros(ndof)
Pbar = np.zeros(ndof)
Pbar[dim*22+1] = -Pref/2
Pbar[dim*23+1] = -Pref/2
# Pbar[dim*22] = Pref/20
# Pbar[dim*23+1] = Pref/20

P0 = P0[free_dof]
Pbar = Pbar[free_dof]

# Start at u = zero, gamma = 0.2
u_n = np.zeros(ndof)
gamma_n = 0.0
gamma_guess = 1

# Initialization Lists for plotting
zeroU = np.zeros(ndof)
u_guess = np.zeros(ndof)
Results = [{'u': u_n, 'gamma': gamma_n, 'R_list': 0, 'Iterations': 0}]
u_list = np.array(u_n)
Uu = [0]
Uv = [0]
G = [0]

# Initial Step to find delta_s^2
# Set Displacements and Construct Global K Matrix
K_global = np.zeros([ndof, ndof])
Fint = np.zeros(ndof)
for ele in T:
    xi = ele.nodeID[0]
    xj = ele.nodeID[1]
    xi_idx = range(xi*dim, (xi+1)*dim)
    xj_idx = range(xj*dim, (xj+1)*dim)
    ele.setDisp(u_n[xi_idx], u_n[xj_idx])
    ke = ele.get_ke()
    K_global[xi*dim:(xi+1)*dim:1, xi*dim:(xi+1)*dim:1] += ke
    K_global[xj*dim:(xj+1)*dim:1, xj*dim:(xj+1)*dim:1] += ke
    K_global[xi*dim:(xi+1)*dim:1, xj*dim:(xj+1)*dim:1] -= ke
    K_global[xj*dim:(xj+1)*dim:1, xi*dim:(xi+1)*dim:1] -= ke
    Fint[xi_idx] -= ele.get_F()
    Fint[xj_idx] += ele.get_F()

Kf = K_global[free_dof, :][:, free_dof]
Fintf = Fint[free_dof]

P = P0 + gamma_guess * Pbar
R = P - Fintf

u_guess[free_dof] = u_n[free_dof] + np.dot(np.linalg.inv(Kf), R)
s2 = np.dot(u_guess[free_dof], u_guess[free_dof]) + a*gamma_guess*gamma_guess
print('(delta_s)^2 = {}, delta_s = {}'.format(s2, np.sqrt(s2)))

flag = 0
for i in range(500):
    for count in range(max_iteration + 1):
        K_global = 0*K_global
        Fint = 0*Fint
        for ele in T:
            xi = ele.nodeID[0]
            xj = ele.nodeID[1]
            xi_idx = range(xi * dim, (xi + 1) * dim)
            xj_idx = range(xj * dim, (xj + 1) * dim)
            ele.setDisp(u_guess[xi_idx], u_guess[xj_idx])
            ke = ele.get_ke()
            K_global[xi * dim:(xi + 1) * dim:1, xi * dim:(xi + 1) * dim:1] += ke
            K_global[xj * dim:(xj + 1) * dim:1, xj * dim:(xj + 1) * dim:1] += ke
            K_global[xi * dim:(xi + 1) * dim:1, xj * dim:(xj + 1) * dim:1] -= ke
            K_global[xj * dim:(xj + 1) * dim:1, xi * dim:(xi + 1) * dim:1] -= ke
            Fint[xi_idx] -= ele.get_F()
            Fint[xj_idx] += ele.get_F()

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

        # Exit if top node hits the ground
        if u_guess[-1] + X[23][1] <= -0.05:
            flag = 1

    if flag == 1:
        break

    # Store Converged Results
    u_list = np.vstack((u_list, u_guess))
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
# Plot of Deformations
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1, figsize=(8, 8))

# Plotting Node Trajectory
for i in range(24):
    plt.plot(u_list[2*i] + X[i][0], u_list[2*i+1] + X[i][1], label="Node {}".format(i), linestyle='--')

# Plotting Final Deformed Shape
for ele in T:
    xi = ele.nodeID[0]
    xj = ele.nodeID[1]
    xi_idx = range(xi * dim, (xi + 1) * dim)
    xj_idx = range(xj * dim, (xj + 1) * dim)
    qi = u_list[xi_idx, -1] + X[xi]
    qj = u_list[xj_idx, -1] + X[xj]
    plt.plot([qi[0], qj[0]], [qi[1], qj[1]], color='black', marker='.')

plt.xlabel('u [m]')
plt.ylabel('v [m]')

plt.axis('equal')
plt.grid(True)
# plt.savefig("CESG506/Code/Assignments/Assignment 3/Problem3-3/Prob3-3_DeformedShape")
plt.show()

# # ----------------------------------------------------------------------------------------------------------------------
# # Plot of gamma vs displacement, part 3 of assignment problem
# # ----------------------------------------------------------------------------------------------------------------------
# plt.figure(1, figsize=(16, 8))
# plt.plot(Uu, G, label="$U_u$", marker='.')
# plt.plot(Uv, G, label="$U_v$", marker='.')
#
# plt.xlabel('displacement [m]')
# plt.ylabel('$\gamma$')
#
# plt.legend(loc='upper center', ncol=1, framealpha=1)
# plt.axhline(y=0, color='black')
# plt.grid(True)
# # plt.savefig("CESG506/Code/Assignments/Assignment 3/Problem3-3/Prob3-3_GammaVsDisp")
# plt.show()

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
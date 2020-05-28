import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, 'CESG506/Code/Assignments/Assignment 6')
from curvedForceBeamClass import *
from helperFunctions import *


# Problem 5-4 Buckling of a beam

# User input
W = 5       # m
H = W/100    # m
EA = 10000  # kN
EI = 10     # kN*m2
Pc = np.pi*np.pi*EI/(W*W)

nEle = 8    # Equally spaced elements

# Convergence criteria and arc-length initiation parameters
max_iteration = 20
tol = 1e-9
a = 0.0
delta_s_Factor = 5
gamma_guess = 0.1

# User enters the H(x) and its derivatives
def elevation(x):
    c = 4*H/W**2
    h = c*x*(W-x)
    dh = c*(2*x-W)
    ddh = c*2
    return (h, dh, ddh)



# ----------------------------------------------------------------------------------------------------
# END USER INPUT (unless you want to change boundary conditions, the following creates a pin-roller beam)
# ----------------------------------------------------------------------------------------------------
dim = 3             # DOF per node
nNds = nEle+1       # Number of nodes
ndof = dim*nNds     # Number of DOFs globally
umidIdx = round(nEle/2)*dim + 1     # Keeping track of
uq1Idx = round(nEle/4)*dim + 1
uq3Idx = round(3*nEle/4)*dim + 1

Lx = W / nEle
Xinit = np.linspace(0, W, nNds)
Yinit = elevation(Xinit)[0]

# Generating Nodes
nodes = []  # [[ID, [x, y], [fixity e.g. [1,1,0] = pin], ...]
fixity = [1, 1, 0]
for idx in range(0, nNds):
    if idx == nNds-1:
        fixity = [0, 1, 0]
    nodes.append([idx, np.array([Xinit[idx], Yinit[idx]]), fixity])
    fixity = [0, 0, 0]

# Generating Elements
elements = []
for i in range(nEle):
    elements.append(curvedForceBeamClass(i, [i, i+1], nodes[i][1], nodes[i+1][1], EA, EI))

# Collect indexes for free and fixed DOFs
free_dof = []
fixed_dof = []
idx = 0
for node in nodes:
    for fixity in node[2]:
        if fixity == 0:
            free_dof.append(idx)
        else:
            fixed_dof.append(idx)
        idx += 1

# Create Load Vector by integrating distributed load on beam shape functions
PrefGlobal = np.zeros(ndof)
PrefGlobal[-3] = -Pc

Pbar = PrefGlobal[free_dof]  # Reference Force Vector based off distributed load


# ----------------------------------------------------------------------------------------------------
# Solve
# ----------------------------------------------------------------------------------------------------
# Initialization Lists for plotting
u_n = np.zeros(ndof)
gamma_n = 0.0
zeroU = np.zeros(ndof)
u_guess = np.zeros(ndof)
Results = [{'q': u_n, 'gamma': gamma_n, 'R_list': 0, 'Iterations': 0}]
u_list = np.array(u_n)
G = [0]
detKt = []

# Initial Step to find delta_s^2
# Set Displacements and Construct Global K Matrix
[K_global, Fint] = beamAssemble(elements, u_n, dim, ndof)
Kf = K_global[free_dof, :][:, free_dof]
detKt.append(np.linalg.det(Kf))

P = gamma_guess * Pbar
R = P - Fint[free_dof]

u_guess[free_dof] = u_n[free_dof] + np.dot(np.linalg.inv(Kf), R)

s2 = np.dot(u_guess[free_dof], u_guess[free_dof]) + a*gamma_guess*gamma_guess
s2 = s2*np.sqrt(delta_s_Factor)

print('(delta_s)^2 = {}, delta_s = {}'.format(s2, np.sqrt(s2)))

flag = 0
for i in range(500):
    for count in range(max_iteration + 1):
        [K_global, Fint] = beamAssemble(elements, u_guess, dim, ndof)

        # Construct Kg Matrix
        Kf = K_global[free_dof, :][:, free_dof]
        Kg = np.hstack((Kf, np.array([-Pbar]).T))
        Kg = np.vstack((Kg, np.append(-2*(u_guess-u_n)[free_dof], -2*a*(gamma_guess-gamma_n))))

        # Construct Rg Vector
        P = gamma_guess * Pbar
        R = P - Fint[free_dof]
        g = np.dot((u_guess-u_n)[free_dof], (u_guess-u_n)[free_dof]) + a*(gamma_guess-gamma_n)*(gamma_guess-gamma_n) - s2
        Rg = np.append([R], [g])

        # Solve for Change in u and gamma guess
        dU = np.dot(np.linalg.inv(Kg), Rg)

        # Update Guess
        u_guess[free_dof] = u_guess[free_dof] + dU[0:len(free_dof)]
        gamma_guess = gamma_guess + dU[-1]

        # Exit if passes Pcr by a factor of 2
        if gamma_guess > 2:
            flag = 1

        if np.linalg.det(Kf) < 0:
            flag = 1
            print('Det(Kf) < 0 at P = {}'.format(gamma_guess*Pc))

        # Check Convergence
        if np.linalg.norm(Rg) <= tol:
            print('u_mid = {}m, P = {:.6f} kN '
                  'Iteration Count = {}, u_{}'.format(u_n[umidIdx], gamma_n * Pc, count, i))
            break

        # Save and Print Results if it does not converge by max iteration
        if count == max_iteration:  # Exit loop if max iteration is reached before convergence
            print('u_mid = {}m, P = {:.6f} kN '
                  'Stopped at Max Iteration = {}'.format(u_guess, gamma_n * Pc, count))
            flag = 1

    if flag == 1:
        break

    # Store Converged Results
    u_list = np.vstack((u_list, u_guess))
    G.append(gamma_guess)
    detKt.append(np.linalg.det(Kf))

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
# Plot of gamma vs det(Kf), part 1 of assignment problem 4
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1, figsize=(8, 8))
plt.plot(detKt, G, label="$det(K_f)$", marker='.')

plt.title('{} Element Mesh, Imperfect Beam'.format(nEle))
plt.xlabel('det(Kt) [m]')
plt.ylabel('$\gamma$')

plt.legend(loc='best', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)
plt.savefig("CESG506/Code/Assignments/Assignment 6/Problem6-4/Prob6-4-2_GammaVsDetKt_{}Ele".format(nEle))
plt.show()


# ----------------------------------------------------------------------------------------------------------------------
# Plot of gamma vs displacement, part 4 of assignment problem 3
# ----------------------------------------------------------------------------------------------------------------------
plt.figure(1, figsize=(8, 8))
plt.plot(u_list[umidIdx], G, label="$U_u$ [y-dir Mid Node]", marker='.')
plt.plot(u_list[uq1Idx], G, label="$U_u$ [y-dir 1/4 Node]", marker='.')
plt.plot(u_list[uq3Idx], G, label="$U_u$ [y-dir 3/4 Node]", marker='.')

plt.title('{} Element Mesh, Imperfect Beam'.format(nEle))
plt.xlabel('displacement [m]')
plt.ylabel('$\gamma$')

plt.legend(loc='best', ncol=1, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)
plt.savefig("CESG506/Code/Assignments/Assignment 6/Problem6-4/Prob6-4-2_GammaVsDisp_{}Ele".format(nEle))
plt.show()




x = 2


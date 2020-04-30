# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 1 - Problem 1-1:
# Study of different formulations for large deformation problems
# Tatsuhiko Sweet, EDITED for assignment 3 on 4/29/2020

import numpy as np
import matplotlib.pyplot as plt

w = 5.5    # m
H = 0.5    # m
EA = 2100  # kN

umax = 2.5*H      # plot for u = [0, umax]
n_pts = 101       # number of points to plot

# Undeformed length parameter calculations
LL = w*w + H*H    # Undeformed length squared
L  = np.sqrt(LL)  # Undeformed length
Ny = H/L          # y-component of Undeformed N vector
k_ln = (EA/L)*Ny*Ny

# Initialize Deformed Parameters
u = np.linspace(0.0, umax, n_pts)    # Displacement Domain
ny = n_pts*[0]    # list of y-component of deformed unit normal vector

# Initialize Strains
# (not necessary to keep for problem statement but is interesting to plot on their own)
e_d = n_pts*[0]

# Initialize P for each case
# Naming convention => P_(strain index, a,b,c, or d)(N = undeformed, n = deformed)
P_dn = n_pts*[0]

# Compute Deformed ny, strains, and P's for each u
for idx in range(n_pts):
    ui = u[idx]
    ll = w*w + (H-ui)*(H-ui)
    l = np.sqrt(ll)
    ny[idx] = (H - ui) / l
    e_d[idx] = 0.5 * np.log(ll / LL)
    P_dn[idx] = -EA * e_d[idx] * ny[idx]

# Plot P's
plt.plot(-u, P_dn, label="$P_{dn}$")
plt.ylabel('P [kN]')
plt.xlabel('u [m]')
plt.legend(loc='best', ncol=3, framealpha=1)
plt.axhline(y=0, color='black')
plt.grid(True)
plt.show()

# Write to text file (deltaU = 0.0125, from 0 to 1.2375)
f = open("CESG506/Code/Assignments/Assignment 3/Problem3-1/FromHW1/HW1outputs", "w+")
for Pi in P_dn:
    f.write(str(Pi))
    f.write('\n')
f.close()


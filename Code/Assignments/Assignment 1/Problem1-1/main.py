# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 1 - Problem 1-1:
# Study of different formulations for large deformation problems
# Tatsuhiko Sweet, 4/14/2020

import numpy as np
import matplotlib.pyplot as plt

w = 5.5    # m
H = 0.5    # m
EA = 2100  # kN

whatToPlot = 'P'  # use 'P' for force, 'e' for strains
umax = 2.5*H      # plot for u = [0, umax]
n_pts = 1001      # number of points to plot

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
e_a = n_pts*[0]
e_b = n_pts*[0]
e_c = n_pts*[0]
e_d = n_pts*[0]

# Initialize P for each case
# Naming convention => P_(strain index, a,b,c, or d)(N = undeformed, n = deformed)
P_ln = n_pts*[0]
P_aN = n_pts*[0]
P_an = n_pts*[0]
P_bN = n_pts*[0]
P_bn = n_pts*[0]
P_cN = n_pts*[0]
P_cn = n_pts*[0]
P_dN = n_pts*[0]
P_dn = n_pts*[0]

# Compute Deformed ny, strains, and P's for each u
for idx in range(n_pts):
    ui = u[idx]
    ll = w*w + (H-ui)*(H-ui)
    l = np.sqrt(ll)
    ny[idx] = (H - ui) / l

    e_a[idx] = (l - L) / np.sqrt(LL)        # Strains
    e_b[idx] = 0.5 * (ll - LL) / LL
    e_c[idx] = 0.5 * (ll - LL) / ll
    e_d[idx] = 0.5 * np.log(ll / LL)

    P_ln[idx] = k_ln*ui                     # Applied force P (just y comp)
    P_aN[idx] = -EA * e_a[idx] * Ny
    P_an[idx] = -EA * e_a[idx] * ny[idx]
    P_bN[idx] = -EA * e_b[idx] * Ny
    P_bn[idx] = -EA * e_b[idx] * ny[idx]
    P_cN[idx] = -EA * e_c[idx] * Ny
    P_cn[idx] = -EA * e_c[idx] * ny[idx]
    P_dN[idx] = -EA * e_d[idx] * Ny
    P_dn[idx] = -EA * e_d[idx] * ny[idx]

# Plotting
# Plot P's
if whatToPlot == 'P':
    plt.plot(u, P_ln, label="$P_{linear}$")
    plt.plot(u, P_an, label="$P_{an}$")
    plt.plot(u, P_bn, label="$P_{bn}$")
    plt.plot(u, P_cn, label="$P_{cn}$")
    plt.plot(u, P_dn, label="$P_{dn}$")
    plt.plot(u, P_aN, label="$P_{aN}$")
    plt.plot(u, P_bN, label="$P_{bN}$")
    plt.plot(u, P_cN, label="$P_{cN}$")
    plt.plot(u, P_dN, label="$P_{dN}$")
    plt.ylabel('P [kN]')
    fileName = 'force_vs_disp'

# Plot Strains
elif whatToPlot == 'e':
    plt.plot(u, e_a, label="$\epsilon_a$")
    plt.plot(u, e_b, label="$\epsilon_b$")
    plt.plot(u, e_c, label="$\epsilon_c$")
    plt.plot(u, e_d, label="$\epsilon_d$")
    plt.ylabel('Strain')
    fileName = 'strain_vs_disp'

plt.xlabel('u [m]')
plt.legend(loc='best', ncol=3, framealpha=1)
plt.xlim([0, umax])
plt.ylim([-1.5, 2])
plt.axhline(y=0, color='black')
plt.grid(True)
plt.savefig("CESG506/Code/Assignments/Assignment 1/Problem1-1/{}".format(fileName))
plt.show()
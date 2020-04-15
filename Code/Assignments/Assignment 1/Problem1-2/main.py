# CESG 506 - Nonlinear Analysis of Structural Systems
# Assignment 1 - Problem 1-2:
# Two degree of freedom (2 DOF) problem - Load Control
# Tatsuhiko Sweet, 4/14/2020

import numpy as np
import matplotlib.pyplot as plt
from helperFunctions import *

# Define geometry and material
w1 = 5.5 # m
w2 = 4.0 # m
H  = 0.5 # m
EA = 2100 # kN

# Convergence criteria
max_iteration = 1000
tol = 1e-12

# P vectors
Pcr = np.array([0, -0.98171345])
gamma = [0, 0.25, 0.5, 0.75, 0.99, 0.999]
# gamma = np.linspace(0.999999, 1.00000001, 2000)  # Used for finding Pcr

# End Inputs ##########################################################################################################

# Undeformed Length Vectors
L1 = np.array([w1, H])
L2 = np.array([-w2, H])

# Results are collect as a list of dictionaries
results = []
#   each dictionary contains:
#   'P' = P vector,
#   'u' = converged displacement vector,
#   'R' = a list of the magnitude of the residual vector at each iteration step up to convergence,
#   'Step' = a list keeping track of which iteration step a residual calculation is on

step = 1
u = np.array([0, 0])
print('Tolerance = {}, s.t. tol >= |R|'.format(tol))

for g in gamma:
    P = g * Pcr
    R_list = []
    step_list = []
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
        R = P - f1 - f2                 # Compute residual vector

        R_list.append(np.linalg.norm(R))
        step_list.append(step)
        step += 1

        if np.linalg.norm(R) <= tol:        # Check for convergence using tol >= |R|
            print('u = {}m, (x,y) = {}m, P = {:.9f} kN = {}Pcr, '
                  'Iteration Count = {}'.format(u, L1 + u, P[1], g, count))
            break
        if count == max_iteration:          # Exit loop if max iteration is reached before convergence
            print('u = {}m, (x,y) = {}m, P = {} kN = {}Pcr, '
                  'Stopped at Max Iteration = {}'.format(u, L1 + u, P[1], g, count))
        u = u + np.dot(np.linalg.inv(Kf), R)    # Update displacement vector
    results.append({'P': P, 'u': u, 'R': R_list, 'Step': step_list, 'Gamma': g})

# Create a table summarizing residuals for each iteration step (Part 4 of Assignment) #################################
#   each row begins with: gamma*Pcr, Py, residuals ...
try:
    f = open("CESG506/Code/Assignments/Assignment 1/Problem1-2/residuals.txt", "w")
    for result in results:
        f.write(str(result['Gamma']) + 'Pcr, ')
        f.write(str(result['P'][1]) + ', ')
        f.write(str(result['R']).strip('[').strip(']') + '\n')
    f.write('Note: Tolerance = {}, Pcr = {} kN'.format(tol, Pcr))
    f.close()
except:
    print("Couldn't write Residuals to a .txt file")

# Plotting Error Results (Part 5 of Assignment) #######################################################################
plot2 = plt.figure(2)
for result in results:
    plt.plot(result['Step'], result['R'], '-o',
             label="$\gamma = ${}".format(result['Gamma']))
plt.plot([], ' ', label='Tol = {}'.format(tol))
plt.axhline(y=tol, color='black')
plt.yscale('log')
plt.ylabel('|R| [kN]')
plt.xlabel('Cummulative Step Count')
plt.legend(loc='upper left', ncol=1, bbox_to_anchor=(1, 1))
plt.grid(True)
plt.tight_layout()
plt.savefig("CESG506/Code/Assignments/Assignment 1/Problem1-2/R_vs_IterationStep")
plt.show()
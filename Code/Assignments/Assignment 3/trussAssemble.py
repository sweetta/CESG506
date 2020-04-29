import numpy as np


def trussAssemble(Trusses, u, dim, ndof):
    K_global = np.zeros([ndof, ndof])
    Fint = np.zeros(ndof)
    for ele in Trusses:
        xi = ele.nodeID[0]
        xj = ele.nodeID[1]
        xi_idx = range(xi * dim, (xi + 1) * dim)
        xj_idx = range(xj * dim, (xj + 1) * dim)
        ele.setDisp(u[xi_idx], u[xj_idx])
        ke = ele.get_ke()
        K_global[xi*dim:(xi+1)*dim:1, xi*dim:(xi+1)*dim:1] += ke
        K_global[xj*dim:(xj+1)*dim:1, xj*dim:(xj+1)*dim:1] += ke
        K_global[xi*dim:(xi+1)*dim:1, xj*dim:(xj+1)*dim:1] -= ke
        K_global[xj*dim:(xj+1)*dim:1, xi*dim:(xi+1)*dim:1] -= ke
        Fint[xi_idx] -= ele.get_F()
        Fint[xj_idx] += ele.get_F()
    return [K_global, Fint]
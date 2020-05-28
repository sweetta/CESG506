import numpy as np


def beamAssemble(elements, q, dim, ndof):
    K_global = np.zeros([ndof, ndof])
    Fint = np.zeros(ndof)
    for ele in elements:
        xi = ele.nodeID[0]
        xj = ele.nodeID[1]
        xi_idx = range(xi * dim, (xi + 1) * dim)
        xj_idx = range(xj * dim, (xj + 1) * dim)
        ele.setDisp(q[xi_idx], q[xj_idx])
        ke = ele.get_ke()
        K_global[xi*dim:xi*dim+6:1, xi*dim:xi*dim+6:1] += ke
        # K_global[xi*dim:(xi+1)*dim:1, xi*dim:(xi+1)*dim:1] += ke
        # K_global[xj*dim:(xj+1)*dim:1, xj*dim:(xj+1)*dim:1] += ke
        # K_global[xi*dim:(xi+1)*dim:1, xj*dim:(xj+1)*dim:1] -= ke
        # K_global[xj*dim:(xj+1)*dim:1, xi*dim:(xi+1)*dim:1] -= ke
        Fint[xi_idx] += ele.get_F()[0:dim]
        Fint[xj_idx] += ele.get_F()[3:6]
    return [K_global, Fint]

# Given the number of integration points, return the location x and the weight w
# return [(x1, w1), ... (xn, wn)]
def gaussQuadValues(n):
    if n == 1:
        return [(0.0, 2.0)]

    if n == 2:
        return [(-0.577350269189626, 1.0),
                 (0.577350269189626, 1.0)]

    if n == 3:
        return [(-0.774596669241483, 0.555555555555556),
                 (0.0, 0.888888888888889),
                 (0.774596669241483, 0.555555555555556)]

    if n == 4:
        return [(-0.339981043584856, 0.652145154862546),
                (-0.861136311594053, 0.347854845137454),
                 (0.861136311594053, 0.347854845137454),
                 (0.339981043584856, 0.652145154862546)]

    if n == 5:
        return [(-0.906179845938664, 0.236926885056189),
                (-0.538469310105683, 0.478628670499366),
                 (0, 0.568888888888889),
                (0.538469310105683, 0.478628670499366),
                 (0.906179845938664, 0.236926885056189)]

def gaussQuad(n, Le):
    xgw = gaussQuadValues(n)    # interval -1 to 1
    xlw = []
    for xw in xgw:
        xi = 0.5*Le*(xw[0]+1)
        wi = 0.5*Le*xw[1]
        xlw.append((xi, wi))
    return xlw

# Shape Functions and their derivatives
def Nu1(x, Le):
    return 1 - x/Le

def Nu2(x, Le):
    return x/Le

def Nv1(x, Le):
    a = x/Le
    return 1-3*a**2+2*a**3

def Nv2(x, Le):
    a = x/Le
    return Le*(a-2*a**2+a**3)

def Nv3(x, Le):
    a = x/Le
    return 3*a**2 - 2*a**3

def Nv4(x, Le):
    a = x/Le
    return Le*(-a**2 + a**3)

def dNu1(x, Le):
    return -1/Le

def dNu2(x, Le):
    return 1/Le

def dNv1(x, Le):
    a = x/Le
    return (6*a**2 - 6*a)/Le

def dNv2(x, Le):
    a = x/Le
    return 3*a**2 - 4*a + 1

def dNv3(x, Le):
    a = x/Le
    return (-6*a**2 + 6*a)/Le

def dNv4(x, Le):
    a = x/Le
    return 3*a**2 - 2*a

def ddNv1(x, Le):
    return 12*x/Le**3 - 6/Le**2

def ddNv2(x, Le):
    return 6*x/Le**2 - 4/Le

def ddNv3(x, Le):
    return 6/Le**2 - 12*x/Le**3

def ddNv4(x, Le):
    return -2/Le + 6*x/Le**2

def dNu(x, Le):
    return np.array([dNu1(x, Le), 0, 0, dNu2(x, Le), 0, 0])

def Nv(x, Le):
    return np.array([0, Nv1(x, Le), Nv2(x, Le), 0, Nv3(x, Le), Nv4(x, Le)])

def dNv(x, Le):
    return np.array([0, dNv1(x, Le), dNv2(x, Le), 0, dNv3(x, Le), dNv4(x, Le)])

def ddNv(x, Le):
    return np.array([ddNv1(x, Le), ddNv2(x, Le), ddNv3(x, Le), ddNv4(x, Le)])

def Beps(x, Le, q, dh):
    dv = np.dot(dNv(x, Le), q)
    return np.array([dNu1(x, Le), (dh+dv)*dNv1(x, Le), (dh+dv)*dNv2(x, Le),
                     dNu2(x, Le), (dh+dv)*dNv3(x, Le), (dh+dv)*dNv4(x, Le)])

def Bphi(x, Le):
    return np.array([0, ddNv1(x, Le), ddNv2(x, Le),
                     0, ddNv3(x, Le), ddNv4(x, Le)])









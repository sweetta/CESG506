import sys
import numpy as np
from helperFunctions import *

class curvedBeamClass(object):
    """
    class: curvedBeamClass

    variables:
        self.ID
        self.X1 = coordinate of node i in undeformed state
        self.X2 = coordinate of node j in undeformed state
        self.EA = EA
        self.EI = EI
        self.Le = undeformed length
        self.Lx = undeformed length in x direction only
        self.dh = linearly interpolated approximation of h'(x)

        self.U1 = displacement of node i (ui, vi, theta_i)
        self.U2 = displacement of node j (uj, vj, theta_j)
        self.strain = axial strain (Henkey)
        self.Fint = internal force vector (np.array)
        self.ke = element stiffness matrix (ndof x ndof np.array)

    methods:
        setDisp(self, U1, U2) ... Recompute deformed state variables from given displacements
                                  at nodes i (U1, np.array) and node j (U2, np.array)
        setFint(self) ........... Internal method to set Force vector (np.array)
        setStiffness(self) ...... Internal method to set ke (6x6 np.array)
        get_Fint(self) .......... External method to return self.F (np.array)
        get_ke(self) ............ External method to return self.ke (6x6 np.array)
    """
    def __init__(self, eleNum, nodeID, X1, X2, EA, EI):
        self.eleNum = eleNum
        self.nodeID = nodeID
        self.X1 = X1
        self.X2 = X2
        self.EA = EA
        self.EI = EI

        self.Le = np.linalg.norm(X2-X1)
        self.Lx = abs(X2[0] - X1[0])
        self.dh = (X2[1] - X1[1])/self.Lx

    def setDisp(self, U1, U2):
        self.ke = np.zeros((6, 6))
        self.Fint = np.zeros(6)
        self.U1 = U1
        self.U2 = U2
        self.q = np.concatenate((U1, U2))
        XW = gaussQuad(3, self.Lx)
        for xw in XW:
            xg = xw[0] + self.X1[0]     # Global x for node i
            Be = Beps(xw[0], self.Lx, self.q, self.dh)
            Bp = Bphi(xw[0], self.Lx)
            du0 = np.dot(dNu(xw[0], self.Lx), self.q)
            dv0 = np.dot(dNv(xw[0], self.Lx), self.q)
            F = self.EA*(du0 + dv0*self.dh + 0.5*dv0*dv0)
            M = self.EI*np.dot(Bp, self.q)
            ka = self.EA*np.outer(Be, Be)
            kb = self.EI*np.outer(Bp, Bp)
            kc = F*np.outer(dNv(xw[0], self.Lx), dNv(xw[0], self.Lx))
            self.ke += xw[1]*(ka+kb+kc)
            self.Fint += xw[1]*(F*Be + M*Bp)


    def get_F(self):
        return self.Fint

    def get_ke(self):
        return self.ke

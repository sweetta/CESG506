import sys
import numpy as np

class trussClass(object):
    """
    class: trussClass

    variables:
        self.ndof = number of DOF, aka dimension of problem
        self.X1 = coordinate of node i in undeformed state
        self.X2 = coordinate of node j in undeformed state
        self.EA = EA
        self.L = undeformed L vector
        self.L2 = undeformed length squared (scalar)

        self.U1 = displacement of node i
        self.U2 = displacement of node j
        self.strain = axial strain (Henkey)
        self.l = deformed l vector
        self.le = deformed length (scalar)
        self.n = normal vector in deformed state
        self.f = internal axial force (scalar)
        self.F = internal axial force in the deformed n direction (np.array)
        self.ke = element stiffness matrix (ndof x ndof np.array)

    methods:
        setDisp(self, U1, U2) ... Recompute deformed state variables from given displacements
                                  at nodes i (U1, np.array) and node j (U1, np.array)
        setStrain(self) ......... Internal method to set strain
        setStiffness(self) ...... Internal method to set ke
        set_n(self) ............. Internal method to set n
        get_F(self) ............. External method to return self.F (np.array)
        get_f(self) ............. External method to return self.f (scalar)
        get_ke(self) ............ External method to return self.ke (ndof x ndof np.array)
    """
    def __init__(self, nodeID, X1, X2, EA):

        self.ndof = len(X1)
        self.X1 = X1
        self.X2 = X2
        self.EA = EA
        self.L = self.X2-self.X1
        self.L2 = np.dot(self.L, self.L)
        self.nodeID = nodeID

        self.U1 = []
        self.U2 = []
        self.strain = 0
        self.l = self.L
        self.le = np.linalg.norm(self.l)
        self.n = self.L/self.le
        self.f = 0
        self.F = 0*self.n
        self.setStiffness()

    def setDisp(self, U1, U2):
        self.U1 = U1
        self.U2 = U2
        self.l = self.L + U2 - U1
        self.le = np.linalg.norm(self.l)
        self.setStrain()
        self.set_n()
        self.f = self.EA*self.strain
        self.F = self.f*self.n
        self.setStiffness()

    def setStrain(self):
        stretchSquared = np.dot(self.l, self.l)/self.L2
        self.strain = 0.5*np.log(stretchSquared)

    def setStiffness(self):
        dim = self.ndof
        f_int = self.f
        n = self.n
        length = self.le
        ke = (self.EA / length) * np.outer(n, n) + \
             (f_int / length) * (np.identity(dim) - np.outer(n, n))
        self.ke = ke

    def set_n(self):
        self.n = self.l/self.le

    def get_F(self):
        return self.F

    def get_f(self):
        return self.f

    def get_ke(self):
        return self.ke






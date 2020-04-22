import sys
import numpy as np

class trussClass(object):
    def __init__(self, nodeID, X1, X2, EA):

        self.ndof = len(X1)
        self.X1 = X1
        self.X2 = X2
        self.EA = EA
        self.L = self.X2-self.X1
        self.L2 = np.dot(self.L, self.L)

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






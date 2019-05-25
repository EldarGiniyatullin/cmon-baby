import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
from scipy import linalg as ln
from numpy import random as rd
from numpy.random import binomial


def signification(x):
    return 2*x-1


def signGen(size):
    return np.apply_along_axis(signification, 0, binomial(1, 0.5, size))


def signGenMatrix(_n, _m):
    B = []
    for i in range(_n):
        B.append(signGen(_m))
    return np.asarray(B)


class array4SPS:
    def __init__(self, _parameters, _realTheta, generator, *args):
        self.parameters = _parameters.copy()
        self.y = []
        for i in range(len(self.parameters)):
            self.y.append(np.dot(self.parameters[i], _realTheta) + generator(args))
        # self.dim = (_parameters[0].shape if _dim is None else _dim)
        self.dim = _parameters[0].shape[0]
        print("Parameters dimension =", self.dim)


# class SPSCalculationResults:
#   def __init__(self):


class SPSmodule:
    def __init__(self, _array4me, _q, _m, _n=None, _begin=0):
        self.y = _array4me.y.copy()
        self.parameters = _array4me.parameters.copy()
        self.q = _q
        self.m = _m
        self.begin = _begin
        self.dim = _array4me.dim
        #TODO: test this one!
        self.n = _array4me.parameters.size if _n is None else _n
        self.oneNth = 1/self.n
        #TODO: test this one!
        self.A = np.transpose(np.append([np.ones(self.n)], signGenMatrix(self.m, self.n), axis = 0))
        self.Rn = np.zeros([self.dim, self.dim])
        self.Rn_root = np.zeros([self.dim, self.dim])
        self.Rn_inv_root = np.zeros([self.dim, self.dim])
        #TODO: test this one!
        self.H = np.full((self.m, self.dim), np.NINF)
        self.S = np.full((self.m, self.dim), np.NINF)
        self.Hnorm = np.full(self.m, np.NINF)
        self.Snorm = np.full(self.m, np.NINF)
        self.g = np.full((self.n, self.dim), np.NINF)
        self.log = []

    def generateA(self):
        self.A = np.transpose(np.append([np.ones(self.n)], signGenMatrix(self.m, self.n), axis = 0))
        return

    def generateAfrom1(self):
        self.A = np.delete(np.append(self.A, [self.A[1]], axis=0), 1, 0)

    def countRn(self):
        self.Rn = np.zeros([self.dim, self.dim])
        for i in range(self.begin, self.begin + self.n):
            self.Rn += np.outer(self.parameters[i], self.parameters[i])
        self.Rn_root = ln.cholesky(self.Rn)
        try:
            self.Rn_inv_root = ln.inv(self.Rn_root)
        except ln.LinAlgError:
            print("Rn is not invertible")
        except ln.ValueError:
            print("Rn is not a matrix!")
        else:
            # cause of cholesky returns L*L.H,
            # but we need the 2nd matrix
            return

    def refreshParameters(self, _parameters):
        self.parameters = _parameters.copy()
        self.countRn()
        return

    def countG(self, t, theta):
        self.g[t] = np.transpose(self.parameters[self.begin + t])*(self.y[self.begin + t] - self.parameters[self.begin + t]*np.transpose(theta))
        return

    def countAllGs(self, theta):
        for t in range(self.n):
            self.countG(t, theta)
        return

    def countH(self, i, theta):
        SPSum = np.zeros([self.dim])
        for t in range(self.begin, self.begin + self.n):
            SPSum += self.A[t][i]*self.g[t]
        return SPSum

    #TODO: not needed in isThetaInRegion()
    def countS(self, i, theta):
        #TODO: precount g_t(theta) (count once and then use saved results)
        return (1/self.n)*(self.Rn_inv_root)*(self.countH(i, theta))

    def isThetaInRegion(self, theta):
        counterS = 0
        counterH = 0
        equalsS = 0
        equalsH = 0
        self.countAllGs(theta)
        self.generateA()
        for i in range(self.m):
            self.H[i] = self.countH(i, theta)
            self.S[i] = self.oneNth * np.matmul(self.Rn_inv_root, (np.transpose(self.H[i])))
            self.Hnorm[i] = np.dot(self.H[i], self.H[i])
            self.Snorm[i] = np.dot(self.S[i], self.S[i])
            if (self.Hnorm[i] < self.Hnorm[0]):
                counterH += 1
            if (self.Snorm[i] < self.Snorm[0]):
                counterS += 1
            if (self.Hnorm[i] == self.Hnorm[0]):
                equalsH += 1
            if (self.Snorm[i] == self.Snorm[0]):
                equalsS += 1
        return True if (counterS + (equalsS // 2)) < (self.m - self.q) else False

    def isThetaInRegionH(self, theta):
        counterH = 0
        equalsH = 0
        for i in range(self.m):
            self.H[i] = self.countH(i, theta)
            self.Hnorm[i] = np.dot(self.H[i], np.transpose(self.H[i]))
            if (self.Hnorm[i] < self.Hnorm[0]):
                counterH += 1
            if (self.Hnorm[i] == self.Hnorm[0]):
                equalsH += 1
        return True if (counterH + (equalsH // 2)) < (self.m - self.q) else False

    def isThetaInRegionHPrev1(self, theta):
        counterH = 0
        equalsH = 0

        # DOES NOT NEED RECOUNT OF A, BEGIN
        self.generateAfrom1()
        self.begin = self.begin + 1
        nplus1 = self.begin + self.n
        gnplus1 = np.transpose(self.parameters[nplus1])*(self.y[nplus1] - self.parameters[nplus1]*np.transpose(theta))
        deltaG = gnplus1 - self.g[0]
        self.g = np.delete(np.append(self.g, [gnplus1], axis=0), 0, 0)
        for i in range(self.m):
            self.H[i] = self.H[i] + deltaG
            self.Hnorm[i] = self.Hnorm[i] + 2*self.A[self.n-1][i]*np.dot(self.H[0], deltaG) # + else!!
            if (self.Hnorm[i] < self.Hnorm[0]):
                counterH += 1
            if (self.Hnorm[i] == self.Hnorm[0]):
                equalsH += 1
        return True if (counterH + (equalsH // 2)) < (self.m - self.q) else False
# def determiner of theta's p-region "in-being"??
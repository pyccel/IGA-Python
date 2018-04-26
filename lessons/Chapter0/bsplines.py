# -*- coding: UTF-8 -*-
#! /usr/bin/python
import numpy as np
from scipy.interpolate import splev
import matplotlib.pyplot as plt

def make_knots(n_elements,p,kind="open"):
    """
    creates a knot vector
    """
    if kind == "open":
        T = np.linspace(0.,1.,n_elements+1)
        T = [0.] *p  + list(T) + [1.] * p
        T = np.array(T)
    else:
        raise NotImplemented("not implemented")

    return T

class Bspline(object):
    def __init__(self, T, p):
        """
        initialize splines for given knot sequence t and degree p
        """
        self.T = T
        self.p = p
        self.N = len(T) - p - 1
        self.c = np.zeros(self.N)

    def __call__(self, x, i=None, n_deriv=0):
        """
        evaluate b-spline starting at node i at x
        """
        if i is not None:
            c = np.zeros_like(self.T)
            if i < 0:
                c[self.N + i]=1.
            else:
                c[i]=1.
        else:
            c = self.c

        tck = (self.T, c, self.p)
        return splev(x, tck, der=n_deriv)

    def greville(self):
        """
        Returns the Greville points
        """
        p = self.p
        T = self.T

        # TODO implement a pure python function and not use igakit
        from igakit import igalib
        return igalib.bsp.Greville(p, T)

    def plot(self, nx=100):
        """
        Plots all splines constructed from a knot sequence
        """
        T = self.T
        p = self.p
        N = self.N
        x = np.linspace(0.0,1.0,nx)

        y = np.zeros((N,nx), dtype=np.double)
        for i in range(0,N):
            y[i]=self(x, i=i)
            plt.plot(x,y[i])

class UniformBspline(Bspline):
    # TODO add xmin and xmax as arguments for the interval bouondaries
    def __init__(self, p, N):
        """
        creates a uniform bspline
        """
        L = list(range(-p, N + p + 1))
        T = np.array(L, dtype=np.float)

        Bspline.__init__(self, T, p)

class CardinalBspline(UniformBspline):
    def __init__(self, p):
        """
        creates a uniform bspline
        """
        N = p + 1
        UniformBspline.__init__(self, p, N)

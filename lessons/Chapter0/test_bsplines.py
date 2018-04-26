# -*- coding: UTF-8 -*-
#! /usr/bin/python

from bsplines import Bspline
from bsplines import CardinalBspline
from bsplines import UniformBspline
import numpy as np
import matplotlib.pyplot as plt


def test_openknot_bspline():
    n = 5
    p = 2
    T = np.linspace(0.,1.,n+1)
    T = [0.] *p  + list(T) + [1.] * p
    T = np.array(T)
    bsp = Bspline(T,p)

    nx = 100
    N = len(T) - p - 1
    x = np.linspace(0.0,1.0,nx)

    y = np.zeros((N,nx), dtype=np.double)
    for i in range(0,N):
        y[i]=bsp(x, i=i)
        plt.plot(x,y[i])
    plt.show()

def test_cardinal_bspline():
    p = 3
    bsp = CardinalBspline(p)

    nx = 200
    x = np.linspace(0, p+1, nx)

    y = np.zeros(nx, dtype=np.double)
    y=bsp(x, i=3)
    plt.plot(x,y)
    plt.show()

def test_cardinal_bspline2():
    p = 3
    for p in [1,2,3,4]:
        bsp = CardinalBspline(p)

        nx = 200
        x = np.linspace(0, p+1, nx)

        y = np.zeros(nx, dtype=np.double)
        y=bsp(x, i=p)
        plt.plot(x,y, label="$p="+str(p)+"$")
    plt.title("Cardinal B-Splines for degrees 1,2,3 and 4.")
    plt.legend()
    plt.show()

def test_uniform_bspline():
    p = 3
    N = 7
    bsp = UniformBspline(p, N)

    nx = 200
    x = np.linspace(0.0,N*1.0,nx)

    y = np.zeros((N,nx), dtype=np.double)
    for i in range(0,N):
        y[i]=bsp(x, i=i)
        plt.plot(x,y[i])
    plt.show()

def test_uniform_bspline2():
    p = 3
    N = 6
    bsp = UniformBspline(p, N)

    xmin = 0. ; xmax = N

    data = [[0,p-1,"--"], [p-1,p,"-"], [p,N,"--"]]
    for d in data:
        xmin = d[0] ; xmax = d[1] ; c = d[2]
        nx = 100
        x = np.linspace(xmin,xmax,nx)
        y = np.zeros((N,nx), dtype=np.double)
        for i in range(0,N):
            y[i]=bsp(x, i=i)
            plt.plot(x,y[i],c+"k")
    plt.title("Non-vanishing Cardinal B-Splines on the interval $[2,3]$")
    plt.legend()
    plt.show()

def test_uniform_bspline3():
    p = 3
    N = 12
    bsp = UniformBspline(p, N)

    nx = 100
    c = np.random.random(N)**4
    plt.ylim(c.min(), c.max())
    for i in range(0,p):
        xmin = i-1
        xmax = i

        x = np.linspace(xmin,xmax,nx)
        y = np.zeros(nx, dtype=np.double)
        for i in range(0,N):
            y += c[i] * bsp(x, i=i)
        plt.plot(x,y)
    plt.title("A Cubic Cardinal B-Spline serie.")
    plt.legend()
    plt.show()

def greville(T, p):
    N = len(T) - p - 1
    t = np.zeros(N)
    for i in range(0, N):
        t[i] = np.sum(T[i:i+p+1])
    t = t / p
    return t

def test_uniform_bspline_10():
    p = 2
    N = 7
    bsp = UniformBspline(p, N)

    t = greville(bsp.T, p)
    c = t**2
    bsp.c = c

    nx = 50
    x = np.linspace(0.0,N*1.0,nx)

    y=bsp(x)
    plt.plot(x,y)
    plt.show()

##########################################
if __name__ == "__main__":
    test_cardinal_bspline()
    test_cardinal_bspline2()
    test_uniform_bspline()
    test_uniform_bspline2()
    test_uniform_bspline3()
    test_openknot_bspline()

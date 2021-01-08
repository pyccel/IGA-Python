from numpy import empty
import numpy as np
from matplotlib import pyplot as plt

# ==========================================================
def find_span( knots, degree, x ):
    # Knot index at left/right boundary
    low  = degree
    high = 0
    high = len(knots)-1-degree

    # Check if point is exactly on left/right boundary, or outside domain
    if x <= knots[low ]: returnVal = low
    elif x >= knots[high]: returnVal = high-1
    else:
        # Perform binary search
        span = (low+high)//2
        while x < knots[span] or x >= knots[span+1]:
            if x < knots[span]:
                high = span
            else:
                low  = span
            span = (low+high)//2
        returnVal = span

    return returnVal

# ==========================================================
def all_bsplines( knots, degree, x, span ):
    left   = empty( degree  , dtype=float )
    right  = empty( degree  , dtype=float )
    values = empty( degree+1, dtype=float )

    values[0] = 1.0
    for j in range(0,degree):
        left [j] = x - knots[span-j]
        right[j] = knots[span+1+j] - x
        saved    = 0.0
        for r in range(0,j+1):
            temp      = values[r] / (right[r] + left[j-r])
            values[r] = saved + right[r] * temp
            saved     = left[j-r] * temp
        values[j+1] = saved

    return values

# ==========================================================
def point_on_bspline_curve(knots, P, x):
    degree = len(knots) - len(P) - 1
    d = P.shape[-1]

    span = find_span( knots, degree, x )
    b    = all_bsplines( knots, degree, x, span )

    c = np.zeros(d)
    for k in range(0, degree+1):
        c[:] += b[k]*P[span-degree+k,:]
    return c

# ==========================================================
def plot_field_1d(knots, degree, u, nx=101, color='b'):
    n = len(knots) - degree - 1

    xmin = knots[degree]
    xmax = knots[-degree-1]

    xs = np.linspace(xmin, xmax, nx)

    P = np.zeros((len(u), 1))
    P[:,0] = u[:]
    Q = np.zeros((nx, 1))
    for i,x in enumerate(xs):
        Q[i,:] = point_on_bspline_curve(knots, P, x)

    plt.plot(xs, Q[:,0], '-'+color)

# ==========================================================
def point_on_bspline_surface(Tu, Tv, P, u, v):
    pu = len(Tu) - P.shape[0] - 1
    pv = len(Tv) - P.shape[1] - 1
    d = P.shape[-1]

    span_u = find_span( Tu, pu, u )
    span_v = find_span( Tv, pv, v )

    bu   = all_bsplines( Tu, pu, u, span_u )
    bv   = all_bsplines( Tv, pv, v, span_v )

    c = np.zeros(d)
    for ku in range(0, pu+1):
        for kv in range(0, pv+1):
            c[:] += bu[ku]*bv[kv]*P[span_u-pu+ku, span_v-pv+kv,:]

    return c

# ==========================================================
def plot_field_2d(knots, degrees, u, nx=101, ny=101):
    T1,T2 = knots
    p1,p2 = degrees

    n1 = len(T1) - p1 - 1
    n2 = len(T2) - p2 - 1

    xmin = T1[p1]
    xmax = T1[-p1-1]

    ymin = T2[p2]
    ymax = T2[-p2-1]

    xs = np.linspace(xmin, xmax, nx)
    ys = np.linspace(ymin, ymax, ny)

    n1,n2 = u.shape

    P = np.zeros((n1, n2, 1))
    P[:,:,0] = u[:,:]
    Q = np.zeros((nx, ny, 1))
    for i1,x in enumerate(xs):
        for i2,y in enumerate(ys):
            Q[i1,i2,:] = point_on_bspline_surface(T1, T2, P, x, y)
    X,Y = np.meshgrid(xs,ys)
    plt.contourf(X, Y, Q[:,:,0])

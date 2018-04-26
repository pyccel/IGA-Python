# coding: utf-8
import numpy as np
from spl.linalg.stencil import StencilVector, StencilMatrix

def assemble_matrix_1d(V, kernel, args=None, M=None):

    # ... sizes
    [s1] = V.vector_space.starts
    [e1] = V.vector_space.ends
    [p1] = V.vector_space.pads

    k1 = V.quad_order
    spans_1 = V.spans
    basis_1 = V.basis
    weights_1 = V.weights
    # ...

    # ... data structure if not given
    if M is None:
        M = StencilMatrix(V.vector_space, V.vector_space)
    # ...

    # ... element matrix
    mat = np.zeros((p1+1,2*p1+1), order='F')
    # ...

    if args is None:
        _kernel = kernel
    else:
        _kernel = lambda p1, k1, bs, w, mat: kernel(p1, k1, bs, w, mat, *args)

    # ... build matrices
    for ie1 in range(s1, e1+1-p1):
        i_span_1 = spans_1[ie1]

        bs = basis_1[:, :, :, ie1]
        w = weights_1[:, ie1]
        _kernel(p1, k1, bs, w, mat)
        s1 = i_span_1 - p1 - 1
        M._data[s1:s1+p1+1,:] += mat[:,:]
    # ...

    return M


def assemble_rhs_1d(V, f):

    # ... sizes
    [s1] = V.vector_space.starts
    [e1] = V.vector_space.ends
    [p1] = V.vector_space.pads

    k1 = V.quad_order
    spans_1 = V.spans
    basis_1 = V.basis
    points_1 = V.points
    weights_1 = V.weights
    # ...

    # ... data structure
    rhs = StencilVector( V.vector_space )
    # ...

    # ... build rhs
    for ie1 in range(s1, e1+1-p1):
        i_span_1 = spans_1[ie1]
        for il_1 in range(0, p1+1):
            i1 = i_span_1 - p1  - 1 + il_1

            v = 0.0
            for g1 in range(0, k1):
                bi_0 = basis_1[il_1, 0, g1, ie1]
                bi_x = basis_1[il_1, 1, g1, ie1]

                x1    = points_1[g1, ie1]
                wvol  = weights_1[g1, ie1]

                v += bi_0 * f(x1) * wvol


            rhs[i1] += v
    # ...

    # ...
    return rhs
    # ...

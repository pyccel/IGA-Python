__all__ = ['assemble_stiffness_2d',
           'assemble_vector_2d'
]

#==============================================================================
def assemble_stiffness_2d(nelements, degree, spans, basis, weights, points, matrix):
    """
    assembling the stiffness matrix using stencil forms
    """

    # ... sizes
    ne1,ne2              = nelements
    p1,p2                = degree
    spans_1, spans_2     = spans
    basis_1, basis_2     = basis
    weights_1, weights_2 = weights
    points_1, points_2   = points

    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
    # ...

    # ... build matrices
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, ne2):
            i_span_2 = spans_2[ie2]
            # evaluation dependant uniquement de l'element

            for il_1 in range(0, p1+1):
                for il_2 in range(0, p2+1):
                    for jl_1 in range(0, p1+1):
                        for jl_2 in range(0, p2+1):
                            i1 = i_span_1 - p1 + il_1
                            j1 = i_span_1 - p1 + jl_1

                            i2 = i_span_2 - p2 + il_2
                            j2 = i_span_2 - p2 + jl_2

                            v = 0.0
                            for g1 in range(0, k1):
                                for g2 in range(0, k2):
                                    bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2]
                                    bi_x = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2]
                                    bi_y = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2]

                                    bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2]
                                    bj_x = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2]
                                    bj_y = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2]

                                    wvol = weights_1[ie1, g1] * weights_2[ie2, g2]

                                    v += (bi_x * bj_x + bi_y * bj_y) * wvol

                            matrix[i1, i2, j1-i1, j2-i2]  += v
    # ...

    return matrix

#==============================================================================
def assemble_vector_2d(f, nelements, degree, spans, basis, weights, points, rhs):
    """
    Assembly procedure for the rhs
    """

    # ... sizes
    ne1,ne2              = nelements
    p1,p2                = degree
    spans_1, spans_2     = spans
    basis_1, basis_2     = basis
    weights_1, weights_2 = weights
    points_1, points_2   = points

    k1 = weights_1.shape[1]
    k2 = weights_2.shape[1]
    # ...

    # ... build rhs
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, ne2):
            i_span_2 = spans_2[ie2]

            for il_1 in range(0, p1+1):
                for il_2 in range(0, p2+1):
                    i1 = i_span_1 - p1 + il_1
                    i2 = i_span_2 - p2 + il_2

                    v = 0.0
                    for g1 in range(0, k1):
                        for g2 in range(0, k2):
                            bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2]
                            bi_x = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2]
                            bi_y = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2]

                            x1    = points_1[ie1, g1]
                            x2    = points_2[ie2, g2]
                            wvol  = weights_1[ie1, g1]*weights_2[ie2, g2]

                            v += bi_0 * f(x1,x2) * wvol

                    rhs[i1,i2] += v
    # ...

    # ...
    return rhs
    # ...

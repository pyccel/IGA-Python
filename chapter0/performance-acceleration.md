# Performance and Acceleration

In this section, we shall see how to accelerate our Python code assembly and get native speed.
We will be using [Numba](https://numba.pydata.org/) and [Pyccel](https://github.com/pyccel/pyccel).

The **StencilMatrix** format are based on a negative indexing, which we provide through a syntactic-sugar approach by overiding the **__getitem__** and **__setitem__** methods. The entries are in fact stored in the the private attribute **\_data** which is a *numpy.NdArray*. We therefor use a shift by adding **p** as described in the **StencilMatrix** method

```python
class StencilMatrix( object ):
    ...
    def _shift_index( index, shift ):
        if isinstance( index, slice ):
            start = None if index.start is None else index.start + shift
            stop  = None if index.stop  is None else index.stop  + shift
            return slice(start, stop, index.step)
        else:
            return index + shift
```

In order to use a Python accelerator, we will first need to pass the *NdArray* as argument to the assembly function and not the **StencilMatrix** object.

Since this will be *ugly* it is highly recommanded that one creates an *interface* function that will be calling the *core* assembly function. It is the later one that we will accelerate using Numba or Pyccel.


```python
from psydac.fem.splines import SplineSpace
from psydac.fem.tensor  import TensorFemSpace
from psydac.linalg.stencil import StencilMatrix
from psydac.linalg.stencil import StencilVector
from psydac.ddm.cart import DomainDecomposition
```

## 1D Case

### Pure Python

The novelties here are

* we add Python Typing Syntax to our function
* we use the matrix as **inout** argument, in the spirit of Fortran. This means that we shall have a procedure and not a function.
* we shift the access to memory when storing the Matrix entries


```python
def assemble_stiffness_1d_pure(nk1: 'int', 
                               p1: 'int',
                               nq1: 'int', 
                               spans_1: 'int[:]', 
                               basis_1: 'double[:,:,:,:]', 
                               weights_1: 'double[:,:]', 
                               points_1: 'double[:,:]', 
                               matrix: 'double[:,:]'):
    """
    assembling the stiffness matrix using stencil forms
    """

    # ... build matrices
    for ie1 in range(0, nk1):
        i_span_1 = spans_1[ie1]
        for il_1 in range(0, p1+1):
            for jl_1 in range(0, p1+1):
                i1 = i_span_1 - p1 + il_1
                j1 = i_span_1 - p1 + jl_1

                v = 0.0
                for g1 in range(0, nq1):
                    bi_0 = basis_1[ie1, il_1, 0, g1]
                    bi_x = basis_1[ie1, il_1, 1, g1]

                    bj_0 = basis_1[ie1, jl_1, 0, g1]
                    bj_x = basis_1[ie1, jl_1, 1, g1]

                    wvol = weights_1[ie1, g1]

                    v += (bi_x * bj_x) * wvol

                # we shift the test index by p1
                matrix[i1, p1+j1-i1]  += v
    # ...

    # NOTE: we will not return the matrix. 
    #       explainations will come later
    #return matrix
```

#### Timing using pure Python


```python
import numpy as np

xmin = 0. ; xmax = 1.
nelements = 400
degree = 3
grid = np.linspace( xmin, xmax, num=nelements+1 )
V = SplineSpace(degree=degree, grid=grid)
dd = DomainDecomposition(ncells=[V.ncells], periods=[False])
V = TensorFemSpace(dd, V)

M = StencilMatrix(V.vector_space, V.vector_space)

# Sizes
[s1] = V.vector_space.starts
[e1] = V.vector_space.ends
[p1] = V.vector_space.pads

# Quadrature data
nk1       = V.quad_grids()[0].num_elements
nq1       = V.quad_grids()[0].num_quad_pts
spans_1   = V.quad_grids()[0].spans
basis_1   = V.quad_grids()[0].basis
points_1  = V.quad_grids()[0].points
weights_1 = V.quad_grids()[0].weights

%timeit assemble_stiffness_1d_pure( nk1, p1, nq1, spans_1, basis_1, weights_1, points_1, M._data )
```

    13.2 ms ± 126 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)


### Using Pyccel


```python
from pyccel.epyccel import epyccel
```


```python
assemble_stiffness_1d_pyccel = epyccel(assemble_stiffness_1d_pure)
```

#### Timing


```python
import numpy as np

xmin = 0. ; xmax = 1.
nelements = 400
degree = 3
grid = np.linspace( xmin, xmax, num=nelements+1 )
V = SplineSpace(degree=degree, grid=grid)
dd = DomainDecomposition(ncells=[V.ncells], periods=[False])
V = TensorFemSpace(dd, V)

M = StencilMatrix(V.vector_space, V.vector_space)

# Sizes
[s1] = V.vector_space.starts
[e1] = V.vector_space.ends
[p1] = V.vector_space.pads

# Quadrature data
nk1       = V.quad_grids()[0].num_elements
nq1       = V.quad_grids()[0].num_quad_pts
spans_1   = V.quad_grids()[0].spans
basis_1   = V.quad_grids()[0].basis
points_1  = V.quad_grids()[0].points
weights_1 = V.quad_grids()[0].weights

%timeit assemble_stiffness_1d_pyccel( nk1, p1, nq1, spans_1, basis_1, weights_1, points_1, M._data )
```

    12.2 µs ± 139 ns per loop (mean ± std. dev. of 7 runs, 100,000 loops each)


## 2D case

### Pure Python


```python
def assemble_stiffness_2d_pure(nk1: 'int', nk2: 'int', 
                               p1:  'int', p2:  'int', 
                               nq1:  'int', nq2:  'int', 
                               spans_1:          'int[:]', spans_2:          'int[:]', 
                               basis_1: 'double[:,:,:,:]', basis_2: 'double[:,:,:,:]', 
                               weights_1:   'double[:,:]', weights_2:   'double[:,:]', 
                               points_1:    'double[:,:]', points_2:    'double[:,:]', 
                               matrix: 'double[:,:,:,:]'):
    """
    assembling the stiffness matrix using stencil forms
    """

    # ... build matrices
    for ie1 in range(0, nk1):
        i_span_1 = spans_1[ie1]
        for ie2 in range(0, nk2):
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
                            for g1 in range(0, nq1):
                                for g2 in range(0, nq2):
                                    bi_0 = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 0, g2]
                                    bi_x = basis_1[ie1, il_1, 1, g1] * basis_2[ie2, il_2, 0, g2]
                                    bi_y = basis_1[ie1, il_1, 0, g1] * basis_2[ie2, il_2, 1, g2]

                                    bj_0 = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 0, g2]
                                    bj_x = basis_1[ie1, jl_1, 1, g1] * basis_2[ie2, jl_2, 0, g2]
                                    bj_y = basis_1[ie1, jl_1, 0, g1] * basis_2[ie2, jl_2, 1, g2]

                                    wvol = weights_1[ie1, g1] * weights_2[ie2, g2]

                                    v += (bi_x * bj_x + bi_y * bj_y) * wvol

                            matrix[i1, i2, p1+j1-i1, p2+j2-i2]  += v
    # ...
```

#### Timing


```python
# create the spline space for each direction
x1min = 0. ; x1max = 1.
nelements1 = 32
degree1 = 3
grid1 = np.linspace( x1min, x1max, num=nelements1+1 )
V1 = SplineSpace(degree=degree1, grid=grid1)

x2min = 0. ; x2max = 1.
nelements2 = 32
degree2 = 3
grid2 = np.linspace( x2min, x2max, num=nelements2+1 )
V2 = SplineSpace(degree=degree2, grid=grid2)

# create the tensor space
dd = DomainDecomposition(ncells=[V1.ncells, V2.ncells], periods=[False, False])
V = TensorFemSpace(dd, V1, V2)

M = StencilMatrix(V.vector_space, V.vector_space)

# Sizes
[s1, s2] = V.vector_space.starts
[e1, e2] = V.vector_space.ends
[p1, p2] = V.vector_space.pads

# Quadrature data
[      nk1,       nk2] = [g.num_elements for g in V.quad_grids()]
[      nq1,       nq2] = [g.num_quad_pts for g in V.quad_grids()]
[  spans_1,   spans_2] = [g.spans        for g in V.quad_grids()]
[  basis_1,   basis_2] = [g.basis        for g in V.quad_grids()]
[ points_1,  points_2] = [g.points       for g in V.quad_grids()]
[weights_1, weights_2] = [g.weights      for g in V.quad_grids()]

%timeit assemble_stiffness_2d_pure( nk1, nk2, p1, p2, nq1, nq2, spans_1, spans_2, basis_1, basis_2, weights_1, weights_2, points_1, points_2, M._data )
```

    5.97 s ± 491 ms per loop (mean ± std. dev. of 7 runs, 1 loop each)


### Using Pyccel


```python
assemble_stiffness_2d_pyccel = epyccel(assemble_stiffness_2d_pure)
```

#### Timing


```python
# create the spline space for each direction
x1min = 0. ; x1max = 1.
nelements1 = 32
degree1 = 3
grid1 = np.linspace( x1min, x1max, num=nelements1+1 )
V1 = SplineSpace(degree=degree1, grid=grid1)

x2min = 0. ; x2max = 1.
nelements2 = 32
degree2 = 3
grid2 = np.linspace( x2min, x2max, num=nelements2+1 )
V2 = SplineSpace(degree=degree2, grid=grid2)

# create the tensor space
dd = DomainDecomposition(ncells=[V1.ncells, V2.ncells], periods=[False, False])
V = TensorFemSpace(dd, V1, V2)

M = StencilMatrix(V.vector_space, V.vector_space)

# Sizes
[s1, s2] = V.vector_space.starts
[e1, e2] = V.vector_space.ends
[p1, p2] = V.vector_space.pads

# Quadrature data
[      nk1,       nk2] = [g.num_elements for g in V.quad_grids()]
[      nq1,       nq2] = [g.num_quad_pts for g in V.quad_grids()]
[  spans_1,   spans_2] = [g.spans        for g in V.quad_grids()]
[  basis_1,   basis_2] = [g.basis        for g in V.quad_grids()]
[ points_1,  points_2] = [g.points       for g in V.quad_grids()]
[weights_1, weights_2] = [g.weights      for g in V.quad_grids()]

%timeit assemble_stiffness_2d_pyccel( nk1, nk2, p1, p2, nq1, nq2, spans_1, spans_2, basis_1, basis_2, weights_1, weights_2, points_1, points_2, M._data )
```

    3.56 ms ± 127 µs per loop (mean ± std. dev. of 7 runs, 100 loops each)


## Exercises

1. Write the accelerated version of the rhs assembly produce using Pyccel and Numba and compute their timing.
2. Perform a benchmark between Pyccel and Numba while varying the B-Spline degree.

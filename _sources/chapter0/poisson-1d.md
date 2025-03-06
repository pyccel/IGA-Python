# B-splines FEM solver for Poisson equation (1D)


Following the previous [section](http://nbviewer.jupyter.org/github/ratnania/IGA-Python/blob/main/lessons/Chapter1/01_introduction_fem.ipynb), we implement here a B-Splines FEM for the Poisson problem in 1D, with homogeneous boundary conditions.

Appropriate data-structures for matrices and vectors, involved when assembling a bilinear and linear form, will be discussed in the next section.

We consider the following Poisson problem in 1D:
$$
-u^{\prime \prime} = f, \quad (0,1)
\\
u(0) = 0  ~~ \mbox{and} ~~ u(1) = 0
$$

Using the algorithms developped in section 1., for assembly a bilinear form and a linear form, we see that we need to assemble

* The Stiffness matrix
* The right hand side (rhs)


```python
# needed imports
from numpy import zeros, ones, linspace, zeros_like, asarray
from matplotlib.pyplot import plot, show
%matplotlib inline
```


```python
# ... assembling the stiffness matrix using stencil forms
def assemble_stiffness(nelements, degree, spans, basis, weights, points, matrix):

    # ... sizes
    ne1       = nelements
    p1        = degree
    spans_1   = spans
    basis_1   = basis
    weights_1 = weights
    points_1  = points
    
    k1 = weights.shape[1]
    # ...

    # ... build matrices
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for il_1 in range(0, p1+1):
            for jl_1 in range(0, p1+1):
                i1 = i_span_1 - p1 + il_1
                j1 = i_span_1 - p1 + jl_1

                v = 0.0
                for g1 in range(0, k1):
                    bi_0 = basis_1[ie1, il_1, 0, g1]
                    bi_x = basis_1[ie1, il_1, 1, g1]                    

                    bj_0 = basis_1[ie1, jl_1, 0, g1]
                    bj_x = basis_1[ie1, jl_1, 1, g1]                    

                    wvol = weights_1[ie1, g1]

                    v += (bi_x * bj_x) * wvol

                matrix[i1, j1]  += v
    # ...

    return matrix    
# ...
```


```python
# ... Assembly procedure for the rhs
def assemble_vector(f, nelements, degree, spans, basis, weights, points, rhs):

    # ... sizes
    ne1       = nelements
    p1        = degree
    spans_1   = spans
    basis_1   = basis
    weights_1 = weights
    points_1  = points
    
    k1 = weights.shape[1]
    # ...

    # ... build rhs
    for ie1 in range(0, ne1):
        i_span_1 = spans_1[ie1]
        for il_1 in range(0, p1+1):
            i1 = i_span_1 - p1 + il_1

            v = 0.0
            for g1 in range(0, k1):
                bi_0 = basis_1[ie1, il_1, 0, g1]
                bi_x = basis_1[ie1, il_1, 1, g1]  

                x1    = points_1[ie1, g1]
                wvol  = weights_1[ie1, g1]

                v += bi_0 * f(x1) * wvol

            rhs[i1] += v
    # ...

    # ...
    return rhs
    # ...
# ...
```

---
**Note**: 
We added the function **f** as an argument of the *assemble_rhs* function. In general, we will avoid calling a function over the quadrature points for performance reasons. 

---

We will need the following imports from the **psydac** library.


```python
# nedded imports
from psydac.core.bsplines import elements_spans  # computes the span for each element
from psydac.core.bsplines import make_knots      # create a knot sequence from a grid
from psydac.core.bsplines import quadrature_grid # create a quadrature rule over the whole 1d grid
from psydac.core.bsplines import basis_ders_on_quad_grid # evaluates all bsplines and their derivatives on the quad grid
from psydac.utilities.quadratures import gauss_legendre # computes gauss legendre quadrature points & weights
from utilities.plot import plot_field_1d # plot a solution for 1d problems
```

Now, let's create a B-spline Finite Element Space $V$ over a grid of $32$ elements, and using **cubic** *B-splines*.

A B-Spline space is defined here through a grid, we then call the function **make_knots** to create an open knot sequence. We finally, compute the spans for every element, to avoid doing it for each quadrature point on the current element, in the assembly procedure.


```python
p  = 3    # spline degree
ne = 32   # number of elements

grid  = linspace(0., 1., ne+1)
knots = make_knots(grid, p, periodic=False)
spans = elements_spans(knots, p)  
```


```python
nelements = len(grid) - 1
nbasis    = len(knots) - p - 1

# we need the value a B-Spline and its first derivative
nderiv = 1

# create the gauss-legendre rule, on [-1, 1]
u, w = gauss_legendre( p )

# for each element on the grid, we create a local quadrature grid
points, weights = quadrature_grid( grid, u, w )

# for each element and a quadrature points, 
# we compute the non-vanishing B-Splines
basis = basis_ders_on_quad_grid( knots, p, points, nderiv, normalization=False )
```

Assembling the **Stiffness** matrix is then done using


```python
stiffness = zeros((nbasis, nbasis))
stiffness = assemble_stiffness(nelements, p, spans, basis, weights, points, matrix=stiffness)
```

### Analytical solution

In finite elements, one of the first tests we do is to check polynomial solutions of with a degree that smaller or equal to the one we chose for our discretization.
Here, we take the analytical solution $u(x) = x(1-x)$ which leads to the rhs $f(x) = - \partial_{xx} u(x) = 2$.


```python
f = lambda x: 2.    
rhs = zeros(nbasis)
rhs = assemble_vector(f, nelements, p, spans, basis, weights, points, rhs=rhs)
```

In order to apply the homogeneous Dirichlet boundary conditions, we only have to remove the first and the last splines, since their values at the boundary is $1$. Due to the partition of unity, the remaining splines are $0$ on the boundary.


```python
# apply homogeneous dirichlet boundary conditions
rhs = rhs[1:-1]
stiffness = stiffness[1:-1, 1:-1]
```

We use a **cg** solver from **scipy** to solve the linear system.


```python
from scipy.sparse.linalg import cg
u, info = cg( stiffness, rhs, rtol=1e-6, maxiter=5000 )
```

Now, we need to add the first and last control points, which are equal to **0**, then proceed to the plot.


```python
u = [0.] + list(u) + [0.]
u = asarray(u)
plot_field_1d(knots, p, u, nx=401)
```


    
![png](images/poisson-1d/output_20_0.png)
    


## The Spline Space class

Since we will be calling over and over the block-lines [5] and [6], we created a class called **SplineSpace** that you can import from **psydac** library.

### Matrix profile

Let's check the profile of our matrix.

```python
import matplotlib.pyplot as plt 

plt.spy(stiffness)
```
    
![png](images/poisson-1d/output_23_1.png)


As expected, it is a sparse matrix and it has a typical profile (band one).

### Next step

In fact, we are not done yet!
The problem now, is that the data structure we are using for the matrix (here as a dense matrix) will consume too much memory, and is not taken into account the locality!

Since on each element, there are exactly $p+1$ consecutive non-vanishing *B-splines*, we know that on each element there are at most $2p+1$ non zeros entries (in every direction). This means that we will need to set entries as the following, in $1D$:

```python
M[i1, j1 - i1]  += v_s
```

in $2D$:

```python
M[i1, i2, j1 - i1, j2 - i2]  += v_s
```

and in $3D$:


```python
M[i1, i2, i3, j1 - i1, j2 - i2, j3 - i3]  += v_s
```

Unfortunatly, this can not be represented as a **numpy.ndarray**, because of the negative indexing!
Moreover, if our aim is to write a parallel code, other consideration must be taken into account.

More details will be found in the next [section](http://nbviewer.jupyter.org/github/ratnania/IGA-Python/blob/main/lessons/Chapter1/03_data_structure.ipynb).


## Exercises

### Exercise 1.

1. Implement a 1D Poisson solver using non homogeneous Dirichlet boundary conditions.

2. Implement a 2D Poisson solver using homogeneous and non homogeneous Dirichlet boundary conditions.

### Exercise 2.

We consider the Heat equation in 1D

$$
\partial_t u = \nabla^2 u
$$

with homogeneous boundary conditions.

1. Write an explicit FEM solver. What do you notice about the dependency of the CFL with respect to the degree? what about the regularity?

2. Write an implicit FEM solver.

### Exercise 3.

We consider a unidirectional flow in the velocity/pressure formulation

$$
\partial u_x = - \frac{1}{\rho} \partial_x p + \nu \partial_{yy}^2 u_x + g_x
$$

where the unknown $u_x$ is a function of $y$ only. 

We will consider $\partial_x p = -2$ and $g_x = 0.4$. The boundary conditions are given by $u_x(y=0) = u_x(y=h) = 0$.

1. Write an explicit and implicit time schemes and solvers.

2. Find numerical CFL in the case of explicit time scheme.

3. Find the steady state.


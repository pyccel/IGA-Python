# TODO - add docstrings

"""
"""

from bsplines    import elements_spans  # computes the span for each element
from bsplines    import make_knots      # create a knot sequence from a grid
from bsplines    import quadrature_grid # create a quadrature rule over the whole 1d grid
from bsplines    import basis_ders_on_quad_grid # evaluates all bsplines and their derivatives on the quad grid
from quadratures import gauss_legendre
from linalg      import StencilVectorSpace

from numpy import linspace

# =================================================================================================
class SplineSpace(object):
    def __init__(self, degree, nelements=None, grid=None, nderiv=1, periodic=False):

        if (nelements is None) and (grid is None):
            raise ValueError('Either nelements or grid must be provided')

        if grid is None:
            grid  = linspace(0., 1., nelements+1)

        knots = make_knots(grid, degree, periodic=periodic)
        spans = elements_spans(knots, degree)

        nelements = len(grid) - 1
        nbasis    = len(knots) - degree - 1

        # create the gauss-legendre rule, on [-1, 1]
        u, w = gauss_legendre( degree )

        # for each element on the grid, we create a local quadrature grid
        points, weights = quadrature_grid( grid, u, w )

        # for each element and a quadrature points,
        # we compute the non-vanishing B-Splines
        basis = basis_ders_on_quad_grid( knots, degree, points, nderiv )

        self._periodic = periodic
        self._knots = knots
        self._spans = spans
        self._grid = grid
        self._degree = degree
        self._nelements = nelements
        self._nbasis = nbasis
        self._points = points
        self._weights = weights
        self._basis = basis

        self._vector_space = StencilVectorSpace([nbasis], [degree], [periodic])

    @property
    def vector_space(self):
        return self._vector_space

    @property
    def periodic(self):
        return self._periodic

    @property
    def knots(self):
        return self._knots

    @property
    def spans(self):
        return self._spans

    @property
    def grid(self):
        return self._grid

    @property
    def degree(self):
        return self._degree

    @property
    def nelements(self):
        return self._nelements

    @property
    def nbasis(self):
        return self._nbasis

    @property
    def points(self):
        return self._points

    @property
    def weights(self):
        return self._weights

    @property
    def basis(self):
        return self._basis

# =================================================================================================
class TensorSpace(object):
    def __init__( self, *args ):
        """."""
        assert all( isinstance( s, SplineSpace ) for s in args )
        self._spaces = tuple(args)

        nbasis   = [V.nbasis   for V in self.spaces]
        degree   = [V.degree   for V in self.spaces]
        periodic = [V.periodic for V in self.spaces]

        self._vector_space = StencilVectorSpace(nbasis, degree, periodic)

    @property
    def vector_space(self):
        return self._vector_space

    @property
    def spaces(self):
        return self._spaces

    @property
    def knots(self):
        return [V.knots for V in self.spaces]

    @property
    def spans(self):
        return [V.spans for V in self.spaces]

    @property
    def grid(self):
        return [V.grid for V in self.spaces]

    @property
    def degree(self):
        return [V.degree for V in self.spaces]

    @property
    def nelements(self):
        return [V.nelements for V in self.spaces]

    @property
    def nbasis(self):
        return [V.nbasis for V in self.spaces]

    @property
    def points(self):
        return [V.points for V in self.spaces]

    @property
    def weights(self):
        return [V.weights for V in self.spaces]

    @property
    def basis(self):
        return [V.basis for V in self.spaces]

# =================================================================================================
def test_1d():
    V = SplineSpace(degree=3, nelements=16)

    M = StencilMatrix(V.vector_space, V.vector_space)
    u = StencilVector(V.vector_space)

def test_2d():
    V1 = SplineSpace(degree=3, nelements=16)
    V2 = SplineSpace(degree=2, nelements=8)
    V = TensorSpace(V1, V2)

    M = StencilMatrix(V.vector_space, V.vector_space)
    u = StencilVector(V.vector_space)

##################################
if __name__ == '__main__':
    from linalg import StencilVector, StencilMatrix

    test_1d()
    test_2d()

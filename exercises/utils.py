import  h5py
import  numpy               as np
import  matplotlib.pyplot   as plt
from    matplotlib          import cm, colors

from sympde.expr            import EssentialBC
from sympde.topology        import element_of

from psydac.api.essential_bc        import apply_essential_bc
from psydac.linalg.basic            import LinearOperator, Vector
from psydac.linalg.block            import BlockVectorSpace
from psydac.linalg.stencil          import StencilVectorSpace


def get_unit_vector(v, n1, n2, n3, pads1, pads2, pads3):

    v *= 0.0
    if n3 is None:
        assert pads3 is None
        if n2 is None:
            raise NotImplementedError('This get_unit_vector method is only implemented in 2D.')
        else:
            v._data[pads1+n1, pads2+n2] = 1.
    else:
        raise NotImplementedError('This get_unit_vector method is only implemented in 2D.')
    
    return v

def toarray(A):
    """Obtain a numpy array representation of a LinearOperator (which has not implemented toarray())."""
    assert isinstance(A, LinearOperator)

    W = A.codomain

    W_is_block = True if isinstance(W, BlockVectorSpace) else False
    if not W_is_block:
        assert isinstance(W, StencilVectorSpace)

    A_arr = np.zeros(A.shape, dtype="float64")
    w = W.zeros()
    At = A.T

    if W_is_block:
        codomain_blocks = [W[i] for i in range(W.n_blocks)]
    else:
        codomain_blocks = (W, )

    start_index = 0
    for k, Wk in enumerate(codomain_blocks):
        w *= 0.
        v = Wk.zeros()
        if len(Wk.npts) == 2:
            npts1, npts2 = Wk.npts
            pads1, pads2 = Wk.pads
            for n1 in range(npts1):
                for n2 in range(npts2):
                    e_n1_n2 = get_unit_vector(v, n1, n2, None, pads1, pads2, None)
                    if W_is_block:
                        w[k] = e_n1_n2
                        e_n1_n2 = w
                    A_n1_n2 = At @ e_n1_n2
                    A_arr[start_index + n1*npts2+n2, :] = A_n1_n2.toarray()
        else:
            raise NotImplementedError('This toarray method is currently only implemented in 2D.')
        start_index += Wk.dimension

    return A_arr

class H1BoundaryProjector2D(LinearOperator):
    def __init__(self, V0, V0_vs, periodic=[False, False]):

        assert all([isinstance(P, bool) for P in periodic])

        self._domain    = V0_vs
        self._codomain  = V0_vs
        self._space     = V0
        self._periodic  = periodic

        self._BC        = self._get_BC()
        
    def copy(self):
        return H1BoundaryProjector2D(self._space, self._domain, self._periodic)

    def _get_BC(self):

        periodic = self._periodic
        if all([P == True for P in periodic]):
            return None
        
        space   = self._space
        u       = element_of(space, name='u')
        bcs     = [EssentialBC(u, 0, side, position=0) for side in space.domain.boundary]

        bcs_x   = [bcs[0], bcs[1]] if periodic[0] == False else []
        bcs_y   = [bcs[2], bcs[3]] if periodic[1] == False else []

        BC      = bcs_x + bcs_y

        return BC

    @property
    def domain(self):
        return self._domain
    
    @property
    def codomain(self):
        return self._codomain
    
    @property
    def dtype(self):
        return None
    
    def tosparse(self):
        raise NotImplementedError
    
    def toarray(self):
        return toarray(self)
    
    def transpose(self, conjugate=False):
        return self

    def dot(self, v, out=None):
        BC = self._BC
        if out is not None:
            assert isinstance(out, Vector)
            assert out.space is self.codomain
        else:
            out = self.codomain.zeros()

        v.copy(out=out)
        apply_essential_bc(out, *BC)
        return out

def plot(gridsize_x, gridsize_y, title, funs, titles, surface_plot=False):

    x = np.linspace(0, 1, gridsize_x+1)
    y = np.linspace(0, 1, gridsize_y+1)
    xx, yy = np.meshgrid(x, y)
    vals = [[] for _ in funs]
    for i, fun in enumerate(funs):
        for xi, yi in zip(xx, yy):
            vals[i].append([fun(xii, yii) for xii, yii in zip(xi, yi)])
        vals[i] = np.array(vals[i])

    n_plots = len(funs)
    if n_plots > 1:
        assert n_plots == len(titles)
    else:
        if titles is not None:
            print('Warning [plot]: will discard argument titles for a single plot')

    fig = plt.figure(figsize=(2.6+4.8*n_plots, 4.8))
    fig.suptitle(title, fontsize=14)

    for i in range(n_plots):
        vmin = np.min(vals[i])
        vmax = np.max(vals[i])
        cnorm = colors.Normalize(vmin=vmin, vmax=vmax)
        ax = fig.add_subplot(1, n_plots, i+1)
        ax.contourf(xx, yy, vals[i], 50, norm=cnorm, cmap='viridis')
        ax.axis('equal')
        fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap='viridis'), ax=ax,  pad=0.05)
        ax.set_xlabel( r'$x$', rotation='horizontal' )
        ax.set_ylabel( r'$y$', rotation='horizontal' )
        if n_plots > 1:
            ax.set_title ( titles[i] )
    plt.show()

    if surface_plot:
        fig = plt.figure(figsize=(2.6+4.8*n_plots, 4.8))
        fig.suptitle(title+' -- surface', fontsize=14)

        for i in range(n_plots):
            vmin = np.min(vals[i])
            vmax = np.max(vals[i])
            cnorm = colors.Normalize(vmin=vmin, vmax=vmax)
            ax = fig.add_subplot(1, n_plots, i+1, projection='3d')
            ax.plot_surface(xx, yy, vals[i], norm=cnorm, cmap='viridis',
                        linewidth=0, antialiased=False)
            fig.colorbar(cm.ScalarMappable(norm=cnorm, cmap='viridis'), ax=ax,  pad=0.05)
            ax.set_xlabel( r'$x$', rotation='horizontal' )
            ax.set_ylabel( r'$y$', rotation='horizontal' )
            if n_plots > 1:
                ax.set_title ( titles[i] )
        plt.show()

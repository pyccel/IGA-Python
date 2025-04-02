import  numpy               as np
import  matplotlib.pyplot   as plt
from    matplotlib          import cm, colors

from scipy.sparse                   import dia_matrix
from scipy.sparse.linalg            import eigsh

from sympde.expr                    import EssentialBC, BilinearForm, integral
from sympde.topology                import element_of, elements_of, Line, Derham
from sympde.topology.datatype       import H1Space, HcurlSpace, L2Space

from psydac.api.discretization      import discretize
from psydac.api.essential_bc        import apply_essential_bc
from psydac.linalg.basic            import LinearOperator, Vector
from psydac.linalg.block            import BlockVectorSpace, BlockLinearOperator
from psydac.linalg.direct_solvers   import BandedSolver
from psydac.linalg.kron             import KroneckerLinearSolver
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

class HcurlBoundaryProjector2D(LinearOperator):
    def __init__(self, V1, V1_vs, periodic=[False, False]):

        assert all([isinstance(P, bool) for P in periodic])

        self._domain    = V1_vs
        self._codomain  = V1_vs
        self._space     = V1
        self._periodic  = periodic

        self._BC        = self._get_BC()
        
    def copy(self):
        return HcurlBoundaryProjector2D(self._space, self._domain, self._periodic)

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

        bcs_x = [bcs[2], bcs[3]] if periodic[1] == False else []
        bcs_y = [bcs[0], bcs[1]] if periodic[0] == False else []

        BC  = [bcs_x, bcs_y]

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
        for outi, BCi in zip(out, BC):
            apply_essential_bc(outi, *BCi)
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

def to_bnd(A):

    dmat = dia_matrix(A.toarray(), dtype=A.dtype)
    la   = abs(dmat.offsets.min())
    ua   = dmat.offsets.max()
    cmat = dmat.tocsr()

    A_bnd = np.zeros((1+ua+2*la, cmat.shape[1]), A.dtype)

    for i,j in zip(*cmat.nonzero()):
        A_bnd[la+ua+i-j, j] = cmat[i,j]

    return A_bnd, la, ua

def matrix_to_bandsolver(A):
    A.remove_spurious_entries()
    A_bnd, la, ua = to_bnd(A)
    return BandedSolver(ua, la, A_bnd)

def get_M1_block_kron_solver_2D(V1, ncells, degree, periodic):
    """
    Given a 2D DeRham sequenece (V0 = H(grad) --grad--> V1 = H(curl) --curl--> V2 = L2)
    discreticed using ncells, degree and periodic,

        domain = Square('C', bounds1=(0, 1), bounds2=(0, 1))
        derham = Derham(domain)
        domain_h = discretize(domain, ncells=ncells, periodic=periodic, comm=comm)
        derham_h = discretize(derham, domain_h, degree=degree),

    returns the inverse of the mass matrix M1 as a BlockLinearOperator consisting of two KroneckerLinearSolvers on the diagonal.
    """
    # assert 3D
    assert len(ncells) == 2
    assert len(degree) == 2
    assert len(periodic) == 2

    # 1D domain to be discreticed using the respective values of ncells, degree, periodic
    domain_1d = Line('L', bounds=(0,1))
    derham_1d = Derham(domain_1d)

    # storage for the 1D mass matrices
    M0_matrices = []
    M1_matrices = []

    # assembly of the 1D mass matrices
    for (n, p, P) in zip(ncells, degree, periodic):

        domain_1d_h = discretize(domain_1d, ncells=[n], periodic=[P])
        derham_1d_h = discretize(derham_1d, domain_1d_h, degree=[p])

        u_1d_0, v_1d_0 = elements_of(derham_1d.V0, names='u_1d_0, v_1d_0')
        u_1d_1, v_1d_1 = elements_of(derham_1d.V1, names='u_1d_1, v_1d_1')

        a_1d_0 = BilinearForm((u_1d_0, v_1d_0), integral(domain_1d, u_1d_0 * v_1d_0))
        a_1d_1 = BilinearForm((u_1d_1, v_1d_1), integral(domain_1d, u_1d_1 * v_1d_1))

        a_1d_0_h = discretize(a_1d_0, domain_1d_h, (derham_1d_h.V0, derham_1d_h.V0))
        a_1d_1_h = discretize(a_1d_1, domain_1d_h, (derham_1d_h.V1, derham_1d_h.V1))

        M_1d_0 = a_1d_0_h.assemble()
        M_1d_1 = a_1d_1_h.assemble()

        M0_matrices.append(M_1d_0)
        M1_matrices.append(M_1d_1)

    V1_1 = V1[0]
    V1_2 = V1[1]

    B1_mat = [M1_matrices[0], M0_matrices[1]]
    B2_mat = [M0_matrices[0], M1_matrices[1]]

    B1_solvers = [matrix_to_bandsolver(Ai) for Ai in B1_mat]
    B2_solvers = [matrix_to_bandsolver(Ai) for Ai in B2_mat]

    B1_kron_inv = KroneckerLinearSolver(V1_1, V1_1, B1_solvers)
    B2_kron_inv = KroneckerLinearSolver(V1_2, V1_2, B2_solvers)

    M1_block_kron_solver = BlockLinearOperator(V1, V1, ((B1_kron_inv, None), 
                                                        (None, B2_kron_inv)))

    return M1_block_kron_solver

def get_eigenvalues(nb_eigs, sigma, A_m, M_m):
    """
    Compute the eigenvalues of the matrix A close to sigma and right-hand-side M
    Function seen and adapted from >>> psydac_dev/psydac/feec/multipatch/examples/hcurl_eigen_pbms_conga_2d.py <<< 
    (Commit a748a4d8c1569a8765f6688d228f65ea6073c252)

    Parameters
    ----------
    nb_eigs : int
        Number of eigenvalues to compute
    sigma : float
        Value close to which the eigenvalues are computed
    A_m : sparse matrix
        Matrix A
    M_m : sparse matrix
        Matrix M
    """

    print()
    print('-----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  ----- ')
    print(
        'computing {0} eigenvalues (and eigenvectors) close to sigma={1} with scipy.sparse.eigsh...'.format(nb_eigs, sigma))
    mode = 'normal'
    which = 'LM'
    ncv = 4 * nb_eigs
    max_shape_splu = 24000
    if A_m.shape[0] >= max_shape_splu:
        raise ValueError(f'Matrix too large.')
        
    eigenvalues, eigenvectors = eigsh(
        A_m, k=nb_eigs, M=M_m, sigma=sigma, mode=mode, which=which, ncv=ncv)

    print("done: eigenvalues found: " + repr(eigenvalues))
    print('-----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----  -----')
    print()

    return eigenvalues, eigenvectors

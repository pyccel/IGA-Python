# TODO - add docstrings

"""
"""

# coding: utf-8
#
# Copyright 2018 Yaman Güçlü

import numpy as np
from scipy.sparse import coo_matrix

__all__ = ['StencilVectorSpace','StencilVector','StencilMatrix']

#===============================================================================
class StencilVectorSpace( object ):
    """
    Vector space for n-dimensional stencil format. Two different initializations
    are possible:

    - serial  : StencilVectorSpace( npts, pads, periods, dtype=float )

    Parameters
    ----------
    npts : tuple-like (int)
        Number of entries along each direction
        (= global dimensions of vector space).

    pads : tuple-like (int)
        Padding p along each direction (number of diagonals is 2*p+1).

    periods : tuple-like (bool)
        Periodicity along each direction.

    dtype : type
        Type of scalar entries.

    """
    def __init__( self, *args, **kwargs ):

        self._init_serial  ( *args, **kwargs )

    # ...
    def _init_serial( self, npts, pads, periods, dtype=float ):

        assert len(npts) == len(pads) == len(periods)

        # Sequential attributes
        self._starts  = tuple( 0   for n in npts )
        self._ends    = tuple( n-1 for n in npts )
        self._pads    = tuple( pads )
        self._periods = tuple( periods )
        self._dtype   = dtype
        self._ndim    = len( npts )

        # Global dimensions of vector space
        self._npts   = tuple( npts )

    #--------------------------------------
    # Abstract interface
    #--------------------------------------
    @property
    def dimension( self ):
        """ The dimension of a vector space V is the cardinality
            (i.e. the number of vectors) of a basis of V over its base field.
        """
        return np.prod( self._npts )

    # ...
    def zeros( self ):
        """
        Get a copy of the null element of the StencilVectorSpace V.

        Returns
        -------
        null : StencilVector
            A new vector object with all components equal to zero.

        """
        return StencilVector( self )

    # ...
    @property
    def npts( self ):
        return self._npts

    # ...
    @property
    def starts( self ):
        return self._starts

    # ...
    @property
    def ends( self ):
        return self._ends

    # ...
    @property
    def pads( self ):
        return self._pads

    # ...
    @property
    def periods( self ):
        return self._periods

    # ...
    @property
    def dtype( self ):
        return self._dtype

    # ...
    @property
    def ndim( self ):
        return self._ndim

#===============================================================================
class StencilVector( object ):
    """
    Vector in n-dimensional stencil format.

    Parameters
    ----------
    V : psydac.linalg.stencil.StencilVectorSpace
        Space to which the new vector belongs.

    """
    def __init__( self, V ):

        assert isinstance( V, StencilVectorSpace )

        sizes = [e-s+2*p+1 for s,e,p in zip(V.starts, V.ends, V.pads)]
        self._sizes = tuple(sizes)
        self._ndim = len(V.starts)
        self._data  = np.zeros( sizes, dtype=V.dtype )
        self._space = V

    #--------------------------------------
    # Abstract interface
    #--------------------------------------
    @property
    def space( self ):
        return self._space

    #...
    def dot( self, v ):

        assert isinstance( v, StencilVector )
        assert v._space is self._space

        res = self._dot(self._data, v._data , self.pads)

        return res

    #...
    @staticmethod
    def _dot(v1, v2, pads):
        ndim = len(v1.shape)
        index = tuple( slice(p,-p) for p in pads)
        return np.dot(v1[index].flat, v2[index].flat)

    #--------------------------------------
    # Other properties/methods
    #--------------------------------------
    @property
    def starts(self):
        return self._space.starts

    # ...
    @property
    def ends(self):
        return self._space.ends

    # ...
    @property
    def pads(self):
        return self._space.pads

    # ...
    def toarray( self ):
        """ return the local array without the padding"""
        idx = tuple( slice(p,-p) for p in self.pads )
        return self._data[idx].flatten()

    # ...
    def __getitem__(self, key):
        index = self._getindex( key )
        return self._data[index]

    # ...
    def __setitem__(self, key, value):
        index = self._getindex( key )
        self._data[index] = value

    #--------------------------------------
    # Private methods
    #--------------------------------------
    def _getindex( self, key ):

        # TODO: check if we should ignore padding elements

        if not isinstance( key, tuple ):
            key = (key,)
        index = []
        for (i,s,p) in zip(key, self.starts, self.pads):
            if isinstance(i, slice):
                start = None if i.start is None else i.start - s + p
                stop  = None if i.stop  is None else i.stop  - s + p
                l = slice(start, stop, i.step)
            else:
                l = i - s + p
            index.append(l)
        return tuple(index)

#===============================================================================
class StencilMatrix( object ):
    """
    Matrix in n-dimensional stencil format.

    This is a linear operator that maps elements of stencil vector space V to
    elements of stencil vector space W.

    For now we only accept V==W.

    Parameters
    ----------
    V : psydac.linalg.stencil.StencilVectorSpace
        Domain of the new linear operator.

    W : psydac.linalg.stencil.StencilVectorSpace
        Codomain of the new linear operator.

    """
    def __init__( self, V, W ):

        assert isinstance( V, StencilVectorSpace )
        assert isinstance( W, StencilVectorSpace )
        assert W.pads == V.pads

        self._pads     = tuple(V.pads)
        dims           = [e-s+2*p+1 for s,e,p in zip(W.starts, W.ends, W.pads)]
        diags          = [2*p+1 for p in self._pads]
        self._data     = np.zeros( dims+diags, dtype=W.dtype )
        self._domain   = V
        self._codomain = W
        self._ndim     = len( dims )

    #--------------------------------------
    # Abstract interface
    #--------------------------------------
    @property
    def domain( self ):
        return self._domain

    # ...
    @property
    def codomain( self ):
        return self._codomain

    # ...
    def tosparse( self, *, with_pads=False ):

        coo = self._tocoo_no_pads()

        return coo

    #--------------------------------------
    # Other properties/methods
    #--------------------------------------

    # ...
    @property
    def pads( self ):
        return self._pads

    # ...
    def __getitem__(self, key):
        index = self._getindex( key )
        return self._data[index]

    # ...
    def __setitem__(self, key, value):
        index = self._getindex( key )
        self._data[index] = value

    # ...
    def _getindex( self, key ):

        nd = self._ndim
        ii = key[:nd]
        kk = key[nd:]

        index = []

        for i,s,p in zip( ii, self._codomain.starts, self._codomain.pads ):
            x = self._shift_index( i, p-s )
            index.append( x )

        for k,p in zip( kk, self._pads ):
            l = self._shift_index( k, p )
            index.append( l )

        return tuple(index)

    # ...
    @staticmethod
    def _shift_index( index, shift ):
        if isinstance( index, slice ):
            start = None if index.start is None else index.start + shift
            stop  = None if index.stop  is None else index.stop  + shift
            return slice(start, stop, index.step)
        else:
            return index + shift

    #...
    def _tocoo_no_pads( self ):

        # Shortcuts
        nr = self._codomain.npts
        nd = self._ndim
        nc = self._domain.npts
        ss = self._codomain.starts
        pp = self._codomain.pads


        ravel_multi_index = np.ravel_multi_index

        # COO storage
        rows = []
        cols = []
        data = []

        # Range of data owned by local process (no ghost regions)
        local = tuple( [slice(p,-p) for p in pp] + [slice(None)] * nd )

        for (index,value) in np.ndenumerate( self._data[local] ):

            # index = [i1-s1, i2-s2, ..., p1+j1-i1, p2+j2-i2, ...]

            xx = index[:nd]  # x=i-s
            ll = index[nd:]  # l=p+k

            ii = [s+x for s,x in zip(ss,xx)]
            jj = [(i+l-p) % n for (i,l,n,p) in zip(ii,ll,nc,self._pads)]

            I = ravel_multi_index( ii, dims=nr, order='C' )
            J = ravel_multi_index( jj, dims=nc, order='C' )

            rows.append( I )
            cols.append( J )
            data.append( value )

        M = coo_matrix(
                (data,(rows,cols)),
                shape = [np.prod(nr),np.prod(nc)],
                dtype = self._domain.dtype
        )

        M.eliminate_zeros()

        return M

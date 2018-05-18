# coding: utf-8

from scipy.sparse import coo_matrix
from numpy import zeros

def coo_from_blocks(matrices, n_block_rows, n_block_cols):

    # ... compute the global nnz
    nnz = 0
    for i in range(0, n_block_rows):
        for j in range(0, n_block_cols):
            nnz += matrices[i][j].nnz
    # ...

    # ... compute number of rows and cols per block
    n_rows = zeros(n_block_rows, dtype=int)
    n_cols = zeros(n_block_cols, dtype=int)

    for i in range(0, n_block_rows):
        n = 0
        for j in range(0, n_block_cols):
            if not(matrices[i][j] is None):
                n = matrices[i][j].shape[0]
                break
        if n == 0:
            raise ValueError('at least one block must be non empty per row')
        n_rows[i] = n

    for j in range(0, n_block_cols):
        n = 0
        for i in range(0, n_block_rows):
            if not(matrices[i][j] is None):
                n = matrices[i][j].shape[1]
                break
        if n == 0:
            raise ValueError('at least one block must be non empty per col')
        n_cols[j] = n
    # ...

    # ...
    data = zeros(nnz)
    rows = zeros(nnz, dtype=int)
    cols = zeros(nnz, dtype=int)
    # ...

    # ...
    n = 0
    for ir in range(0, n_block_rows):
        for ic in range(0, n_block_cols):
            if not(matrices[ir][ic] is None):
                A = matrices[ir][ic]

                n += A.nnz

                shift_row = 0
                if ir > 0:
                    shift_row = sum(n_rows[:ir])

                shift_col = 0
                if ic > 0:
                    shift_col = sum(n_cols[:ic])

                rows[n-A.nnz:n] = A.row[:] + shift_row
                cols[n-A.nnz:n] = A.col[:] + shift_col
                data[n-A.nnz:n] = A.data
    # ...

    # ...
    nr = n_rows.sum()
    nc = n_cols.sum()

    coo = coo_matrix((data, (rows, cols)), shape=(nr, nc))
    coo.eliminate_zeros()
    # ...

    return coo

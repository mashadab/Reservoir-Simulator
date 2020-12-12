from scipy.sparse import csr_matrix, issparse


def spdiaginv(A):
    """
    Compute the inverse of a sparse diagonal matrix

    Parameters
    ----------
    A : (M,M) ndarray or sparse matrix
        square matrix to be inverted

    Returns
    -------
    Ainv : (M,M) ndarray or sparse matrix
        inverse of `A`

    Notes
    -----
    This computes the sparse inverse of `A`.  If the inverse of `A` is expected
    to be non-sparse, it will likely be faster to convert `A` to dense and use
    scipy.linalg.inv.

    Examples
    --------
    >>> from scipy.sparse import csc_matrix
    >>> from scipy.sparse.linalg import inv
    >>> A = csc_matrix([[1., 0.], [1., 2.]])
    >>> Ainv = inv(A)
    >>> Ainv
    <2x2 sparse matrix of type '<class 'numpy.float64'>'
        with 3 stored elements in Compressed Sparse Column format>
    >>> A.dot(Ainv)
    <2x2 sparse matrix of type '<class 'numpy.float64'>'
        with 2 stored elements in Compressed Sparse Column format>
    >>> A.dot(Ainv).todense()
    matrix([[ 1.,  0.],
            [ 0.,  1.]])

    .. versionadded:: 0.12.0

    """
    #check input
    if not issparse(A):
        raise TypeError('Input must be a sparse matrix')

    n = A.shape[0]       #finding the dimension
    
    Ainv = csr_matrix((n,n))
    
    for i in range(0,n):
        Ainv[i,i] = 1.0/A[i,i]
        
    return Ainv
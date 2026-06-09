"""
Dynamic Number Path Integral core computations.

Key formula:
  1 + W_i = sum_{L=0}^{inf} (T^L)_{ii} = ((I - T)^{-1})_{ii}
"""

from mpmath import mp

def identity_matrix(n=9):
    """Return n x n identity matrix."""
    I = [[mp.mpc(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        I[i][i] = mp.mpc(1)
    return I


def matrix_minus(A, B, n=9):
    """Return A - B for n x n matrices."""
    return [[A[i][j] - B[i][j] for j in range(n)] for i in range(n)]


def matrix_inverse(M, n=9):
    """
    Compute inverse of n x n complex matrix using LU decomposition.
    Returns (M_inv, success_flag).
    """
    # Flatten to list for mpmath's lu_solve
    A_flat = [M[i][j] for i in range(n) for j in range(n)]
    
    try:
        # Use mpmath's built-in matrix inverse
        A_mat = mp.matrix([[M[i][j] for j in range(n)] for i in range(n)])
        A_inv = A_mat ** (-1)
        
        result = [[A_inv[i, j] for j in range(n)] for i in range(n)]
        return result, True
    except Exception:
        return None, False


def compute_W(T, n=9):
    """
    Compute closed-path corrections W_i for all domains.

    Parameters
    ----------
    T : 9x9 transfer matrix
    n : matrix dimension (default 9)

    Returns
    -------
    W : dict {vertex_index: W_i}
        Closed-path corrections for each domain.
    """
    I_mat = identity_matrix(n)
    I_minus_T = matrix_minus(I_mat, T, n)
    
    inv, success = matrix_inverse(I_minus_T, n)
    
    if not success:
        return None
    
    W = {}
    for i in range(n):
        # W_i = ((I-T)^{-1})_{ii} - 1
        W[i] = inv[i][i] - mp.mpc(1)
    
    return W


def compute_effective_correction(W, i):
    """
    Compute |1 + W_i|^2 for domain i.
    This is the effective closed-path correction factor.
    """
    one_plus_W = mp.mpc(1) + W[i]
    return abs(one_plus_W) ** 2

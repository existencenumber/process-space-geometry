"""
Dynamic Number Path Integral core computations.
"""

from mpmath import mp

def compute_W(T):
    """
    Compute closed-path corrections W_i for all domains.
    Returns dict {vertex_index: W_i} or None if inversion fails.
    """
    n = 9
    # 直接用 mp.matrix 构造，不要手动转列表
    T_mat = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            T_mat[i, j] = T[i][j]

    I_mat = mp.eye(n)
    M = I_mat - T_mat

    try:
        inv = M ** (-1)
    except Exception:
        return None

    W = {}
    for i in range(n):
        W[i] = inv[i, i] - mp.mpc(1)
    return W


def compute_effective_correction(W, i):
    """Compute |1 + W_i|^2 for domain i."""
    return abs(mp.mpc(1) + W[i]) ** 2

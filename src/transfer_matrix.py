"""
Construction of the 9x9 transfer matrix T.

T_ij = exp(i * gamma_ij * ln(C_j / C_i))  if edge i->j exists
T_ij = 0                                  otherwise

All non-zero matrix elements have unit modulus.
The dependence on alpha enters through complex powers alpha^{±i/4}, etc.
"""

from mpmath import mp
from src.constants import (
    GAMMA_ANALYTIC, GAMMA_LIMIT, GAMMA_TOP, GAMMA_FREE,
    capacities, VERTEX_ORDER
)

def build_transfer_matrix(alpha):
    """
    Construct the 9x9 transfer matrix T for a given alpha.

    Parameters
    ----------
    alpha : mp.mpf
        Fine-structure constant value.

    Returns
    -------
    T : list of lists (9x9 complex matrix)
        Transfer matrix with 28 non-zero entries.
    """
    C = capacities(alpha)
    sqrt2 = mp.sqrt(2)

    # Initialize 9x9 zero matrix
    T = [[mp.mpc(0) for _ in range(9)] for _ in range(9)]

    # Helper: T_ij = exp(i * gamma * ln(C_j / C_i)) = (C_j/C_i)^{i*gamma}
    def set_edge(i, j, gamma):
        ratio = C[VERTEX_ORDER[j]] / C[VERTEX_ORDER[i]]
        T[i][j] = mp.e ** (mp.mpc(0, 1) * gamma * mp.log(ratio))

    # ---- Edges from A (vertex 0) ----
    set_edge(0, 1, GAMMA_FREE)       # A -> M: exponential map
    set_edge(0, 2, GAMMA_LIMIT)      # A -> I: Riemann limit
    set_edge(0, 3, GAMMA_LIMIT)      # A -> D: difference quotient limit
    set_edge(0, 8, GAMMA_FREE)       # A -> C: identity morphism

    # ---- Edges from M (vertex 1) ----
    set_edge(1, 0, GAMMA_ANALYTIC)   # M -> A: logarithm map
    set_edge(1, 2, GAMMA_ANALYTIC)   # M -> I: logarithmic derivative
    set_edge(1, 4, GAMMA_ANALYTIC)   # M -> S: Mellin transform
    set_edge(1, 8, GAMMA_FREE)       # M -> C: identity morphism'

    # ---- Edges from I (vertex 2) ----
    set_edge(2, 0, GAMMA_LIMIT)      # I -> A: inverse Riemann
    set_edge(2, 1, GAMMA_ANALYTIC)   # I -> M: inverse logarithmic derivative
    set_edge(2, 3, GAMMA_FREE)       # I -> D: FTC
    set_edge(2, 4, GAMMA_ANALYTIC)   # I -> S: Laplace transform
    set_edge(2, 5, GAMMA_LIMIT)      # I -> P: functional limit

    # ---- Edges from D (vertex 3) ----
    set_edge(3, 0, GAMMA_LIMIT)      # D -> A: inverse difference quotient
    set_edge(3, 2, GAMMA_FREE)       # D -> I: inverse FTC
    set_edge(3, 4, GAMMA_ANALYTIC)   # D -> S: Fourier transform

    # ---- Edges from S (vertex 4) ----
    set_edge(4, 1, GAMMA_ANALYTIC)   # S -> M: inverse Mellin
    set_edge(4, 2, GAMMA_ANALYTIC)   # S -> I: inverse Laplace
    set_edge(4, 3, GAMMA_ANALYTIC)   # S -> D: inverse Fourier

    # ---- Edges from P (vertex 5) ----
    set_edge(5, 2, GAMMA_LIMIT)      # P -> I: inverse functional limit
    set_edge(5, 6, GAMMA_TOP)        # P -> B: 2D topology

    # ---- Edges from B (vertex 6) ----
    set_edge(6, 5, GAMMA_TOP)        # B -> P: inverse 2D topology
    set_edge(6, 7, GAMMA_TOP)        # B -> H: braid homotopy

    # ---- Edges from H (vertex 7) ----
    set_edge(7, 6, GAMMA_TOP)        # H -> B: inverse braid homotopy
    set_edge(7, 8, GAMMA_TOP)        # H -> C: morphism categorification

    # ---- Edges from C (vertex 8) ----
    set_edge(8, 0, GAMMA_FREE)       # C -> A: identity morphism
    set_edge(8, 1, GAMMA_FREE)       # C -> M: identity morphism'
    set_edge(8, 7, GAMMA_TOP)        # C -> H: inverse morphism categorification

    return T


def count_nonzero(T):
    """Count non-zero entries in the transfer matrix (should be 28)."""
    count = 0
    for i in range(9):
        for j in range(9):
            if abs(T[i][j]) > mp.mpf('1e-50'):
                count += 1
    return count


def verify_unitarity(T):
    """Verify that all non-zero entries have unit modulus."""
    max_dev = mp.mpf(0)
    for i in range(9):
        for j in range(9):
            if abs(T[i][j]) > mp.mpf('1e-50'):
                dev = abs(abs(T[i][j]) - 1)
                if dev > max_dev:
                    max_dev = dev
    return max_dev

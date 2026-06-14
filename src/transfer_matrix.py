"""
Construction of the 9x9 transfer matrix T with regularization.
"""

from mpmath import mp
from src.constants import (
    GAMMA_ANALYTIC, GAMMA_LIMIT, GAMMA_TOP, GAMMA_FREE,
    capacities, VERTEX_ORDER
)

def build_transfer_matrix(alpha, eps=None):
    """
    Construct the 9x9 transfer matrix T.
    eps: regularization parameter (default = 1/28^2 for vacuum fluctuation).
    """
    if eps is None:
        eps = mp.mpf(1) / (28 ** 2)

    C = capacities(alpha)
    T = [[mp.mpc(0) for _ in range(9)] for _ in range(9)]

    def set_edge(i, j, gamma):
        ratio = C[VERTEX_ORDER[j]] / C[VERTEX_ORDER[i]]
        val = mp.e ** (mp.mpc(0, 1) * gamma * mp.log(ratio))
        if gamma != 0:                     # non-free edge: apply decay
            val *= mp.exp(-eps)
        T[i][j] = val

    # Edges from A (0)
    set_edge(0, 1, GAMMA_FREE)      # A->M: exponential map (free)
    set_edge(0, 2, GAMMA_LIMIT)     # A->I: Riemann limit
    set_edge(0, 3, GAMMA_LIMIT)     # A->D: diff. quotient limit
    set_edge(0, 8, GAMMA_FREE)      # A->C: identity morphism (free)

    # Edges from M (1)
    set_edge(1, 0, GAMMA_ANALYTIC)  # M->A: logarithm
    set_edge(1, 2, GAMMA_ANALYTIC)  # M->I: log-derivative
    set_edge(1, 4, GAMMA_ANALYTIC)  # M->S: Mellin transform
    set_edge(1, 8, GAMMA_FREE)      # M->C: identity morphism (free)

    # Edges from I (2)
    set_edge(2, 0, GAMMA_LIMIT)     # I->A: inverse Riemann
    set_edge(2, 1, GAMMA_ANALYTIC)  # I->M: inverse log-derivative
    set_edge(2, 3, GAMMA_FREE)      # I->D: FTC (free)
    set_edge(2, 4, GAMMA_ANALYTIC)  # I->S: Laplace transform
    set_edge(2, 5, GAMMA_LIMIT)     # I->P: functional limit

    # Edges from D (3)
    set_edge(3, 0, GAMMA_LIMIT)     # D->A: inverse diff. quotient
    set_edge(3, 2, GAMMA_FREE)      # D->I: inverse FTC (free)
    set_edge(3, 4, GAMMA_ANALYTIC)  # D->S: Fourier transform

    # Edges from S (4)
    set_edge(4, 1, GAMMA_ANALYTIC)  # S->M: inverse Mellin
    set_edge(4, 2, GAMMA_ANALYTIC)  # S->I: inverse Laplace
    set_edge(4, 3, GAMMA_ANALYTIC)  # S->D: inverse Fourier

    # Edges from P (5)
    set_edge(5, 2, GAMMA_LIMIT)     # P->I: inverse functional limit
    set_edge(5, 6, GAMMA_TOP)       # P->B: 2D topology

    # Edges from B (6)
    set_edge(6, 5, GAMMA_TOP)       # B->P: inverse 2D topology
    set_edge(6, 7, GAMMA_TOP)       # B->H: braid homotopy

    # Edges from H (7)
    set_edge(7, 6, GAMMA_TOP)       # H->B: inverse braid homotopy
    set_edge(7, 8, GAMMA_TOP)       # H->C: morphism categorification

    # Edges from C (8)
    set_edge(8, 0, GAMMA_FREE)      # C->A: identity morphism (free)
    set_edge(8, 1, GAMMA_FREE)      # C->M: identity morphism (free)
    set_edge(8, 7, GAMMA_TOP)       # C->H: inverse categorification

    return T

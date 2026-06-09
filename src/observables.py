"""
Computation of physical constants from the DNPI transfer matrix.
"""

from mpmath import mp
from src.constants import PI, E, GAMMA, ALPHA0, capacities, VERTEX_ORDER
from src.transfer_matrix import build_transfer_matrix
from src.dnpi import compute_W, compute_effective_correction

# Vertex indices
A_VERTEX = 0   # Addition domain
M_VERTEX = 1   # Multiplication domain
B_VERTEX = 6   # Braid domain


def compute_all_constants(alpha):
    """
    Compute all physical constants given a self-consistent alpha.

    Parameters
    ----------
    alpha : mp.mpf
        Self-consistent fine-structure constant.

    Returns
    -------
    results : dict
        Dictionary of computed constants.
    """
    results = {}
    results['alpha'] = alpha
    results['alpha_inv'] = 1 / alpha
    
    # Build transfer matrix and compute W
    T = build_transfer_matrix(alpha)
    W = compute_W(T)
    C = capacities(alpha)
    
    # ---- Weinberg angle ----
    corr_A = compute_effective_correction(W, A_VERTEX)
    corr_M = compute_effective_correction(W, M_VERTEX)
    
    # sin^2 theta_W = corr_M / (corr_M + (1/3) * alpha * corr_A)
    ratio = (mp.mpf(1)/mp.mpf(3)) * alpha * corr_A
    sin2_thetaW = corr_M / (corr_M + ratio)
    results['sin2_thetaW'] = sin2_thetaW
    
    # ---- Strong coupling ----
    corr_B = compute_effective_correction(W, B_VERTEX)
    
    # alpha_{s,0} = alpha * pi^2 * exp(pi - e) * (1 + (pi-e)/(pi+e))
    alpha_s0 = alpha * (PI**2) * mp.exp(PI - E) * (1 + (PI - E)/(PI + E))
    alpha_s = alpha_s0 / corr_B
    results['alpha_s'] = alpha_s
    
    # ---- Muon-electron mass ratio ----
    phi_r = (3 * PI / 2) * alpha
    
    # Angular phase with torsion correction
    # f_ang = |1 + W_{M->S->M}| / |1 + W_M|
    # Approximate by extracting sub-diagonal block contribution
    inv_T = None
    I_mat = [[mp.mpc(0) for _ in range(9)] for _ in range(9)]
    for i in range(9):
        I_mat[i][i] = mp.mpc(1)
    I_minus_T = [[I_mat[i][j] - T[i][j] for j in range(9)] for i in range(9)]
    
    try:
        A_mat = mp.matrix([[I_minus_T[i][j] for j in range(9)] for i in range(9)])
        A_inv = A_mat ** (-1)
        inv_T = [[A_inv[i, j] for j in range(9)] for i in range(9)]
    except Exception:
        inv_T = None
    
    if inv_T is not None:
        # Extract W_{M->S->M} contribution
        # S vertex = 4, M vertex = 1
        W_M_S_M = inv_T[M_VERTEX][4] * inv_T[4][M_VERTEX]
        f_ang = abs(mp.mpc(1) + W_M_S_M) / abs(mp.mpc(1) + W[M_VERTEX])
    else:
        # Fallback: use torsion estimate
        f_ang = 1 - (GAMMA / (2*PI))**2
    
    phi_ang = mp.mpf(1)/mp.mpf(2) * f_ang
    phi_total = mp.sqrt(phi_r**2 + phi_ang**2)
    
    mass_ratio = (1 - mp.cos(phi_total)) / (1 - mp.cos(phi_r))
    results['mass_ratio_mu_e'] = mass_ratio
    
    return results


def print_results(results):
    """Print computed constants alongside experimental values."""
    print("\n" + "=" * 70)
    print("Process Space Geometry - Computed Physical Constants")
    print("=" * 70)
    print(f"{'Constant':<25} {'Theory':<25} {'Experiment':<25}")
    print("-" * 70)
    
    print(f"{'alpha^(-1)':<25} {results['alpha_inv']:<25.12f} {'137.035999084(21)':<25}")
    print(f"{'sin^2(theta_W)':<25} {results['sin2_thetaW']:<25.6f} {'0.23122(3)':<25}")
    print(f"{'alpha_s(M_Z)':<25} {results['alpha_s']:<25.6f} {'0.1180(9)':<25}")
    print(f"{'m_mu/m_e':<25} {results['mass_ratio_mu_e']:<25.6f} {'206.7682830(46)':<25}")
    print("-" * 70)

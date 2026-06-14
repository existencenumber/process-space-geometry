"""
Computation of physical constants from the DNPI transfer matrix.
"""

from mpmath import mp
from src.constants import PI, E, GAMMA, ALPHA0, capacities
from src.transfer_matrix import build_transfer_matrix
from src.dnpi import compute_W, compute_effective_correction

A_VERTEX = 0
M_VERTEX = 1
B_VERTEX = 6

def compute_all_constants(alpha):
    results = {}
    results['alpha'] = alpha
    results['alpha_inv'] = 1 / alpha

    T = build_transfer_matrix(alpha)
    W = compute_W(T)
    C = capacities(alpha)

    corr_A = compute_effective_correction(W, A_VERTEX)
    corr_M = compute_effective_correction(W, M_VERTEX)

    # --- Weinberg angle (correct formula) ---
    # tan^2 θ_W = (1/3) * α * (corr_M / corr_A)
    tan2 = (mp.mpf(1)/3) * alpha * corr_M / corr_A
    sin2 = tan2 / (1 + tan2)
    results['sin2_thetaW'] = sin2

    # --- Strong coupling ---
    corr_B = compute_effective_correction(W, B_VERTEX)
    alpha_s0 = alpha * (PI**2) * mp.exp(PI - E) * (1 + (PI - E)/(PI + E))
    alpha_s = alpha_s0 / corr_B
    results['alpha_s'] = alpha_s

    # --- Muon-electron mass ratio ---
    phi_r = (3 * PI / 2) * alpha
    f_ang = 1 - (GAMMA / (2*PI))**2   # torsion correction to angular phase
    phi_ang = mp.mpf(1)/2 * f_ang
    phi_total = mp.sqrt(phi_r**2 + phi_ang**2)
    mass_ratio = (1 - mp.cos(phi_total)) / (1 - mp.cos(phi_r))
    results['mass_ratio_mu_e'] = mass_ratio

    return results

def print_results(results):
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

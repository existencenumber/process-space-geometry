"""
Unit tests for physical constant computations.
"""

from mpmath import mp
from src.iteration import solve_alpha
from src.observables import compute_all_constants


def test_alpha_value():
    """Alpha should be close to the CODATA value."""
    alpha, _ = solve_alpha(verbose=False)
    alpha_inv = 1 / alpha
    codata = mp.mpf('137.035999084')
    rel_diff = abs(alpha_inv - codata) / codata
    assert rel_diff < mp.mpf('1e-6'), f"Alpha deviates too much: rel_diff = {rel_diff}"


def test_weinberg_angle():
    """Weinberg angle should be in physical range."""
    alpha, _ = solve_alpha(verbose=False)
    results = compute_all_constants(alpha)
    sin2 = results['sin2_thetaW']
    assert sin2 > mp.mpf('0.2'), f"sin^2 theta_W too small: {sin2}"
    assert sin2 < mp.mpf('0.3'), f"sin^2 theta_W too large: {sin2}"


def test_strong_coupling():
    """Strong coupling should be in physical range."""
    alpha, _ = solve_alpha(verbose=False)
    results = compute_all_constants(alpha)
    alphas = results['alpha_s']
    assert alphas > mp.mpf('0.1'), f"alpha_s too small: {alphas}"
    assert alphas < mp.mpf('0.15'), f"alpha_s too large: {alphas}"


def test_mass_ratio():
    """Mass ratio should be close to experimental value."""
    alpha, _ = solve_alpha(verbose=False)
    results = compute_all_constants(alpha)
    ratio = results['mass_ratio_mu_e']
    exp_val = mp.mpf('206.7682830')
    rel_diff = abs(ratio - exp_val) / exp_val
    assert rel_diff < mp.mpf('1e-3'), f"Mass ratio deviates: rel_diff = {rel_diff}"

"""
Unit tests for alpha iteration.
"""

from mpmath import mp
from src.constants import ALPHA0
from src.iteration import solve_alpha, verify_convergence


def test_convergence():
    """Alpha iteration should converge to a stable value."""
    alpha, history = solve_alpha(verbose=False)
    assert len(history) > 0, "No iteration steps recorded"
    assert history[-1][2] > 137.03, f"Alpha^{-1} too small: {history[-1][2]}"
    assert history[-1][2] < 137.04, f"Alpha^{-1} too large: {history[-1][2]}"


def test_consistency():
    """Self-consistent alpha should reproduce itself in one more step."""
    from src.iteration import iteration_step
    alpha, _ = solve_alpha(verbose=False)
    alpha_new = iteration_step(alpha)
    assert alpha_new is not None, "Iteration step failed"
    rel_diff = abs(alpha_new - alpha) / alpha
    assert rel_diff < mp.mpf('1e-10'), f"Not self-consistent: rel_diff = {rel_diff}"


def test_initial_guess_independence():
    """Different initial guesses should converge to the same value."""
    alpha1, _ = solve_alpha(alpha_init=mp.mpf(1)/mp.mpf(137), verbose=False)
    alpha2, _ = solve_alpha(alpha_init=mp.mpf(1)/mp.mpf(100), verbose=False)
    alpha3, _ = solve_alpha(alpha_init=mp.mpf(1)/mp.mpf(200), verbose=False)
    assert abs(alpha1 - alpha2) / alpha1 < mp.mpf('1e-20')
    assert abs(alpha1 - alpha3) / alpha1 < mp.mpf('1e-20')

"""
Unit tests for transfer matrix construction.
"""

from mpmath import mp
from src.constants import ALPHA0
from src.transfer_matrix import build_transfer_matrix, count_nonzero, verify_unitarity
from src.iteration import solve_alpha


def test_nonzero_count():
    """Transfer matrix should have exactly 28 non-zero entries."""
    alpha, _ = solve_alpha(verbose=False)
    T = build_transfer_matrix(alpha)
    assert count_nonzero(T) == 28, f"Expected 28 non-zero entries, got {count_nonzero(T)}"


def test_unitarity():
    """All non-zero entries should have unit modulus."""
    alpha, _ = solve_alpha(verbose=False)
    T = build_transfer_matrix(alpha)
    dev = verify_unitarity(T)
    assert dev < mp.mpf('1e-50'), f"Unitarity violation: {dev}"


def test_matrix_dimensions():
    """Transfer matrix should be 9x9."""
    alpha, _ = solve_alpha(verbose=False)
    T = build_transfer_matrix(alpha)
    assert len(T) == 9
    for row in T:
        assert len(row) == 9

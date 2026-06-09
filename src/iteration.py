"""
Fixed-point iteration for the fine-structure constant alpha.

Self-consistent equation:
  alpha = alpha_0 / |1 + W_M(alpha)|^2
"""

from mpmath import mp
from src.constants import ALPHA0
from src.transfer_matrix import build_transfer_matrix
from src.dnpi import compute_W, compute_effective_correction

# Vertex indices
M_VERTEX = 1   # Multiplication domain

def iteration_step(alpha):
    """
    Perform one step of the fixed-point iteration.

    Parameters
    ----------
    alpha : mp.mpf
        Current value of alpha.

    Returns
    -------
    alpha_new : mp.mpf
        Updated value of alpha.
    """
    T = build_transfer_matrix(alpha)
    W = compute_W(T)
    
    if W is None:
        return None
    
    correction = compute_effective_correction(W, M_VERTEX)
    alpha_new = ALPHA0 / correction
    
    return alpha_new


def solve_alpha(alpha_init=None, tol=mp.mpf('1e-30'), max_iter=100, 
                verbose=True):
    """
    Solve the self-consistent equation for alpha by fixed-point iteration.

    Parameters
    ----------
    alpha_init : mp.mpf or None
        Initial guess (default: 1/137).
    tol : mp.mpf
        Convergence tolerance.
    max_iter : int
        Maximum number of iterations.
    verbose : bool
        Print iteration progress.

    Returns
    -------
    alpha : mp.mpf
        Self-consistent fine-structure constant.
    history : list of (iteration, alpha, alpha_inv)
        Convergence history.
    """
    if alpha_init is None:
        alpha = mp.mpf(1) / mp.mpf(137)
    else:
        alpha = alpha_init
    
    history = []
    
    for n in range(max_iter):
        alpha_new = iteration_step(alpha)
        
        if alpha_new is None:
            if verbose:
                print(f"Iteration {n}: matrix inversion failed")
            break
        
        alpha_inv = 1 / alpha_new
        history.append((n, float(alpha_new), float(alpha_inv)))
        
        if verbose and n <= 12:
            print(f"Iter {n:3d}: alpha^(-1) = {alpha_inv:.15f}")
        
        if abs(alpha_new - alpha) < tol:
            if verbose:
                print(f"\nConverged after {n+1} iterations.")
                print(f"Final alpha^(-1) = {1/alpha_new:.15f}")
            return alpha_new, history
        
        alpha = alpha_new
    
    if verbose:
        print(f"\nWarning: Did not converge within {max_iter} iterations.")
    return alpha, history


def verify_convergence(alpha, n_extra=10):
    """
    Verify that the solution is stable by performing additional iterations.
    """
    for _ in range(n_extra):
        alpha_new = iteration_step(alpha)
        if alpha_new is None:
            return False, alpha
        alpha = alpha_new
    return True, alpha

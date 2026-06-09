"""
Fundamental constants and parameters for Process Space Geometry.

All computations use 100 decimal digits of precision (mpmath).
The transcendental numbers pi and e are initialized to 100-digit
precision. Euler's constant gamma is derived within the theory
from the torsion dynamic equation.
"""

from mpmath import mp, pi, e

# ============================================================
# High-precision settings
# ============================================================
mp.dps = 100  # 100 decimal digits of working precision

# ============================================================
# Fundamental transcendental numbers
# ============================================================
PI = pi          # mp.pi, 100-digit precision
E  = e           # mp.e,  100-digit precision

# Euler constant gamma derived from torsion dynamic equation:
# gamma = 196 * ln(pi/e) / (49 + ln(pi/e))
TAU = mp.log(PI / E)  # ln(pi/e)
GAMMA = 196 * TAU / (49 + TAU)

# ============================================================
# Tree-level fine-structure constant alpha_0
# alpha_0 = (pi - e)^2 / (pi^2 * sqrt(2*pi)) * (1 + gamma^2/(2*pi)^2)
# ============================================================
def compute_alpha0():
    """Compute the tree-level fine-structure constant alpha_0."""
    num = (PI - E) ** 2
    den = PI**2 * mp.sqrt(2 * PI)
    corr = 1 + GAMMA**2 / (2 * PI)**2
    alpha0 = num / den * corr
    return alpha0

ALPHA0 = compute_alpha0()

# ============================================================
# Characteristic exponents
# ============================================================
GAMMA_ANALYTIC = -1    # Analytic maps (exp, log, Fourier, Laplace, Mellin, log-deriv)
GAMMA_LIMIT    = +1    # Limit maps (Riemann, difference quotient, functional integral)
GAMMA_TOP      = +1    # Topological maps (braid, homotopy, categorification)
GAMMA_FREE     =  0    # Free maps (identity correspondences, FTC)

# ============================================================
# Quantum capacities in the broken phase (functions of alpha)
# ============================================================
def capacities(alpha):
    """Return quantum capacities for all 9 domains given alpha."""
    sqrt_alpha = mp.sqrt(alpha)
    return {
        'A': mp.mpf(1),                    # Addition domain
        'M': 1 / sqrt_alpha,               # Multiplication domain
        'I': alpha ** (-mp.mpf(1)/4),      # Integration domain
        'D': alpha ** (-mp.mpf(1)/4),      # Differentiation domain
        'S': 1 / sqrt_alpha,               # Spectral domain (= C_M)
        'P': alpha ** (-mp.mpf(1)/3),      # Path integral domain
        'B': mp.sqrt(2),                   # Braid domain (topological constant)
        'H': mp.mpf(1),                    # Homotopy domain
        'C': mp.mpf(1),                    # Category domain
    }

# ============================================================
# Vertex ordering
# ============================================================
# A=0, M=1, I=2, D=3, S=4, P=5, B=6, H=7, C=8
VERTEX_ORDER = ['A', 'M', 'I', 'D', 'S', 'P', 'B', 'H', 'C']

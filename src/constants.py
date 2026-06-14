"""
Fundamental constants and parameters for Process Space Geometry.
All computations use 100 decimal digits of precision.
"""

from mpmath import mp, pi, e

mp.dps = 100

PI = pi
E  = e

# Euler constant from torsion dynamic equation
TAU = mp.log(PI / E)
GAMMA = 196 * TAU / (49 + TAU)

# Pure tree-level fine-structure constant (no torsion correction)
# α₀ = (π - e)² / (π² √(2π))
ALPHA0 = (PI - E) ** 2 / (PI**2 * mp.sqrt(2 * PI))

# Characteristic exponents
GAMMA_ANALYTIC = -1    # Analytic maps
GAMMA_LIMIT    = +1    # Limit maps
GAMMA_TOP      = +1    # Topological maps
GAMMA_FREE     =  0    # Free maps (identity, FTC)

def capacities(alpha):
    """Quantum capacities in the broken phase."""
    sqrt_alpha = mp.sqrt(alpha)
    return {
        'A': mp.mpf(1),
        'M': 1 / sqrt_alpha,
        'I': alpha ** (-mp.mpf(1)/4),
        'D': alpha ** (-mp.mpf(1)/4),
        'S': 1 / sqrt_alpha,
        'P': alpha ** (-mp.mpf(1)/3),
        'B': mp.sqrt(2),
        'H': mp.mpf(1),
        'C': mp.mpf(1),
    }

VERTEX_ORDER = ['A', 'M', 'I', 'D', 'S', 'P', 'B', 'H', 'C']

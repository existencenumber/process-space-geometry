#!/usr/bin/env python3
"""
Reproduce all numerical results from the Process Space Geometry papers.

Usage:
    python reproduce_all.py

This script performs the complete computation:
    1. Fixed-point iteration for alpha
    2. Computation of Weinberg angle, strong coupling, mass ratio
    3. Comparison with experimental values
    4. Output of verification report

Requirements:
    pip install mpmath
"""

import sys
import os
from mpmath import mp

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from src.constants import PI, E, GAMMA, TAU, ALPHA0
from src.iteration import solve_alpha
from src.observables import compute_all_constants, print_results


def main():
    print("=" * 70)
    print("Process Space Geometry - Complete Numerical Verification")
    print("=" * 70)
    print(f"Working precision: {mp.dps} decimal digits")
    print(f"pi     = {PI}")
    print(f"e      = {E}")
    print(f"ln(pi/e) = {TAU}")
    print(f"gamma  = {GAMMA}")
    print(f"alpha0^(-1) = {1/ALPHA0:.6f}")
    print("=" * 70)
    
    # Step 1: Solve self-consistent alpha
    print("\n[Step 1] Fixed-point iteration for alpha...")
    alpha, history = solve_alpha(verbose=True)
    
    # Step 2: Compute all constants
    print("\n[Step 2] Computing all physical constants...")
    results = compute_all_constants(alpha)
    
    # Step 3: Print results
    print_results(results)
    
    # Step 4: Verification report
    print("\n[Step 3] Verification Report")
    print("-" * 70)
    
    # Compare with experimental values
    exp_values = {
        'alpha_inv': 137.035999084,
        'sin2_thetaW': 0.23122,
        'alpha_s': 0.1180,
        'mass_ratio_mu_e': 206.7682830,
    }
    
    theory_values = {
        'alpha_inv': float(results['alpha_inv']),
        'sin2_thetaW': float(results['sin2_thetaW']),
        'alpha_s': float(results['alpha_s']),
        'mass_ratio_mu_e': float(results['mass_ratio_mu_e']),
    }
    
    for key in exp_values:
        t_val = theory_values[key]
        e_val = exp_values[key]
        dev = (t_val - e_val) / e_val
        status = "PASS" if abs(dev) < 1e-3 else "CHECK"
        print(f"  {key:<20}: theory={t_val:.10f}, exp={e_val:.10f}, "
              f"dev={dev:.2e} [{status}]")
    
    print("-" * 70)
    print("\nVerification complete.")
    print("All results should match those reported in the paper.")
    print("\nFor questions, contact: xwp499913478@gmail.com")


if __name__ == "__main__":
    main()

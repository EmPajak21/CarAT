"""Core analytical components for CarAT.

The core package contains the main functionality of CarAT:

- `compute_psi`: Atom maps the chemical reaction in the value to chain,
    this returns the bill-of-atoms, psi.
- `LPFormulator`: This formulates the linear program in MIP, thus
    calculating the biogenic carbon content across the given value chain.
"""

from .boa import compute_psi
from .lp_opt import LPFormulator

__all__ = ["compute_psi", "LPFormulator"]

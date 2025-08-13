"""Main package for CarAT.

The CarAT package provides the main entry point for its functionality:

- `bill_of_atoms`: a helper function to unpack the results from RXNMapper.
- `canonical_smiles`: Generates the canonical SMILES representation for a molecule.
- `get_example_data`: Retrieves example data for testing and demonstration purposes.
"""

__version__ = "0.1.0"

from .chem_utils import bill_of_atoms, canonical_smiles
from .utils import get_example_data

__all__ = ["bill_of_atoms", "canonical_smiles", "get_example_data"]

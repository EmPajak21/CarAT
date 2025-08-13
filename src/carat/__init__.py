"""Main package for CarAT.

The CarAT package provides the main entry point for its functionality:

- `bill_of_atoms`: a helper function to unpack the results from RXNMapper.
- `canonical_smiles`: Generates the canonical SMILES representation for a molecule.
- `get_example_data`: Retrieves example data for testing and demonstration purposes.


Copyright (C) 2025  Imperial College London / BASF

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

__version__ = "0.1.0"

from .chem_utils import bill_of_atoms, canonical_smiles
from .utils import get_example_data

__all__ = ["bill_of_atoms", "canonical_smiles", "get_example_data"]

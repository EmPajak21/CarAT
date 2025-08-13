"""Unit tests for the `chem_utils` module.

Verifies that:
- `canonical_smiles` canonicalizes or rejects invalid input.
- `unmap_smiles` strips atom-map annotations.
- `smiles2molecule` and `molecule2smiles` convert between SMILES and RDKit Mol.
- `contains_carbon` correctly detects (or excludes) carbon atoms.
"""

import pandas as pd
from rdkit import Chem

from carat.chem_utils import (
    bill_of_atoms,
    canonical_smiles,
    contains_carbon,
    molecule2smiles,
    smiles2molecule,
    unmap_smiles,
)


def test_canonical_smiles():
    """Ensure `canonical_smiles` returns a valid canonical SMILES or None.

    Checks that:
    - A valid SMILES string is returned unchanged if already canonical.
    - A non-canonical SMILES is converted to its canonical form.
    - Invalid or None inputs yield None.
    """
    # Valid SMILES
    assert canonical_smiles("CCO") == "CCO"
    assert canonical_smiles("OCC") == "CCO"

    # Invalid SMILES
    assert canonical_smiles("invalid") is None

    # None input
    assert canonical_smiles(None) is None


def test_unmap_smiles():
    """Verify that `unmap_smiles` removes atom-map labels correctly.

    Checks that:
    - Mapped SMILES lose their map numbers.
    - Unmapped SMILES remain unchanged.
    - Invalid or None inputs yield None.
    """
    # SMILES with atom-map numbers
    assert unmap_smiles("[CH3:1][CH2:2][OH:3]") == "CCO"

    # SMILES without atom-map numbers
    assert unmap_smiles("CCO") == "CCO"

    # Invalid SMILES
    assert unmap_smiles("invalid") is None

    # None input
    assert unmap_smiles(None) is None


def test_smiles2molecule_and_molecule2smiles():
    """Test round-trip conversion between SMILES and RDKit Mol objects.

    Ensures that:
    - A valid SMILES string produces an RDKit Mol.
    - Converting that Mol back yields the original SMILES.
    - Invalid SMILES or None inputs produce None.
    """
    # Valid SMILES
    smi = "CCO"
    mol = smiles2molecule(smi)
    assert isinstance(mol, Chem.Mol)
    assert molecule2smiles(mol) == smi

    # Invalid SMILES
    assert smiles2molecule("invalid") is None
    assert molecule2smiles(None) is None


def test_contains_carbon():
    """Confirm that `contains_carbon` correctly identifies carbon presence.

    Tests that:
    - SMILES with at least one carbon atom return True.
    - SMILES without carbon (e.g., "Cl") return False.
    - None or NaN inputs return False.
    """
    # SMILES with carbon
    assert contains_carbon("CCO") is True

    # SMILES without carbon
    assert contains_carbon("Cl") is False

    # None or NaN input
    assert contains_carbon(None) is False
    assert contains_carbon(float("nan")) is False


def test_smiles_round_trip_conversion():
    """Test round-trip conversion from SMILES to Mol and back to SMILES.

    Verifies that for a variety of valid SMILES strings, converting to an RDKit Mol
    via `smiles2molecule` and then back to SMILES via `molecule2smiles` yields
    the same canonical SMILES representation.
    """
    examples = [
        "CCO",
        "O=C=O",
        "c1ccccc1",
        "C1CCCCC1",
    ]
    for smi in examples:
        # Parse the SMILES into a molecule
        mol = smiles2molecule(smi)
        assert mol is not None, f"Failed to parse SMILES: {smi}"

        # Convert back to SMILES and canonicalize both
        round_trip = molecule2smiles(mol)
        expected = canonical_smiles(smi)
        assert round_trip == expected, (
            f"Round-trip mismatch for {smi!r}: got {round_trip!r}, "
            f"expected {expected!r}"
        )


def test_bill_of_atoms():
    """Test the `bill_of_atoms` function with a simple mapped reaction.

    Verifies that:
    - Atom mappings are correctly parsed from educts and products.
    - The output DataFrame contains the correct atom counts and mappings.
    """
    # Example reaction: CH4 + Cl2 -> CH3Cl
    mapped_rxn = "[CH4:1].[Cl:2][Cl:3]>>[CH3:1][Cl:2]"

    # Expected output DataFrame
    expected_data = [
        ("CCl", "C", "C", 1),
        ("CCl", "C", "H", 3),
        ("CCl", "ClCl", "Cl", 1),
    ]
    expected_df = pd.DataFrame(
        expected_data, columns=["PROD", "EDUCT", "ATOM", "COUNT"]
    )

    # Call the function
    result_df = bill_of_atoms(mapped_rxn)

    # Ensure the result matches the expected DataFrame
    pd.testing.assert_frame_equal(
        result_df.sort_values(by=["PROD", "EDUCT", "ATOM"]).reset_index(drop=True),
        expected_df.sort_values(by=["PROD", "EDUCT", "ATOM"]).reset_index(drop=True),
    )

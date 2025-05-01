import pytest
from chem_utils import canonical_smiles, unmap_smiles, smiles2molecule, molecule2smiles, contains_carbon
from rdkit import Chem

def test_canonical_smiles():
    # Valid SMILES
    assert canonical_smiles("CCO") == "CCO"
    assert canonical_smiles("OCC") == "CCO"  # Canonical form

    # Invalid SMILES
    assert canonical_smiles("invalid") is None

    # None input
    assert canonical_smiles(None) is None

def test_unmap_smiles():
    # SMILES with atom-map numbers
    assert unmap_smiles("[CH3:1][CH2:2][OH:3]") == "CCO"

    # SMILES without atom-map numbers
    assert unmap_smiles("CCO") == "CCO"

    # Invalid SMILES
    assert unmap_smiles("invalid") is None

    # None input
    assert unmap_smiles(None) is None

def test_smiles2molecule_and_molecule2smiles():
    # Valid SMILES
    smi = "CCO"
    mol = smiles2molecule(smi)
    assert isinstance(mol, Chem.Mol)
    assert molecule2smiles(mol) == smi

    # Invalid SMILES
    assert smiles2molecule("invalid") is None
    assert molecule2smiles(None) is None

def test_contains_carbon():
    # SMILES with carbon
    assert contains_carbon("CCO") is True

    # SMILES without carbon
    assert contains_carbon("Cl") is False

    # None or NaN input
    assert contains_carbon(None) is False
    assert contains_carbon(float('nan')) is False
"""Chemical utilities for CarAT.

Helpers for SMILES canonicalization, atom-map cleanup, carbon detection,
and mapped‐reaction atom counting:

- canonical_smiles(smiles) -> Optional[str]
- unmap_smiles(smi) -> Optional[str]
- smiles2molecule(smi) -> Optional[Mol]
- molecule2smiles(mol) -> Optional[str]
- contains_carbon(smiles) -> bool
- update_smiles(row) -> str
- bill_of_atoms(mapped_rxn) -> pd.DataFrame
"""

import re
from typing import Optional, Union

import pandas as pd
from rdkit import Chem
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol


def canonical_smiles(smiles: Optional[str]) -> Optional[str]:
    """
    Return the RDKit‐canonical SMILES string.

    Parameters
    ----------
    smiles : str or None
        Input SMILES string to be canonicalized.

    Returns
    -------
    str or None
        Canonical SMILES if input is valid; otherwise None.
    """
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol) if mol else None


def unmap_smiles(smi: Union[str, Chem.Mol, None]) -> Optional[str]:
    """
    Remove atom‐map numbers from a SMILES string or Mol object.

    Parameters
    ----------
    smi : str, rdkit.Chem.rdchem.Mol, or None
        SMILES string or Mol for which to strip mapping numbers.

    Returns
    -------
    str or None
        Unmapped SMILES if parsing succeeds; otherwise None.
    """
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else smi
    if mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(mol)
    return None


def smiles2molecule(smi: str) -> Union[Mol, None]:
    """
    Convert a SMILES string to an RDKit Mol object.

    Parameters
    ----------
    smi : str
        SMILES representation of the molecule.

    Returns
    -------
    rdkit.Chem.rdchem.Mol or None
        RDKit Mol object if parsing succeeds; otherwise None.
    """
    return MolFromSmiles(smi) if smi is not None else None


def molecule2smiles(mol: Union[Mol, None]):
    """
    Convert an RDKit Mol object to a SMILES string.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol or None
        RDKit Mol to convert.

    Returns
    -------
    str or None
        SMILES representation if Mol is valid; otherwise None.
    """
    return MolToSmiles(mol) if mol is not None else None


def contains_carbon(smiles):
    """
    Check if a SMILES string contains at least one un‐mapped carbon atom.

    This excludes chlorine ('Cl') and NaN entries.

    Parameters
    ----------
    smiles : str or None
        SMILES string to inspect.

    Returns
    -------
    bool
        True if at least one carbon atom (not part of 'Cl') is present; False otherwise.
    """
    if pd.isna(smiles):
        return False
    else:
        return bool(re.search(r"\bC(?!l)\b|C[^l1]", smiles))


def update_smiles(row):
    """
    Normalize SMILES values in a DataFrame row.

    Converts 'None' or NaN to the literal 'unknown', otherwise canonicalizes.

    Parameters
    ----------
    row : pandas.Series
        DataFrame row containing a "SMILES" column.

    Returns
    -------
    str
        'unknown' if input is missing/invalid; otherwise the canonical SMILES.
    """
    if (row["SMILES"] == "None") or pd.isna(row["SMILES"]):
        return "unknown"
    else:
        return canonical_smiles(row["SMILES"])


def bill_of_atoms(mapped_rxn: str) -> pd.DataFrame:
    """
    Build a count of product atoms by their educt origin.

    Parses a mapped reaction SMILES, tracks atom‐map indices from educts,
    and tallies how many atoms in the product derive from each educt.

    Parameters
    ----------
    mapped_rxn : str
        Reaction SMILES with atom maps, formatted as 'educts>>product'.

    Returns
    -------
    pandas.DataFrame
        Aggregated counts with columns ["PROD", "EDUCT", "ATOM", "COUNT"].
    """
    left, prod_mapped = mapped_rxn.split(">>")
    educt_frags = left.split(".")

    # Where does each map index come from?
    origin = {}
    for frag in educt_frags:
        mol = Chem.MolFromSmiles(frag)
        if not mol:
            continue
        educt = unmap_smiles(frag)
        for a in mol.GetAtoms():
            origin[a.GetAtomMapNum()] = educt

    prod_mol = Chem.MolFromSmiles(prod_mapped)
    product = unmap_smiles(prod_mapped)

    rows = []
    for a in prod_mol.GetAtoms():
        src = origin.get(a.GetAtomMapNum())
        if not src:
            continue  # atom was not present in any educt
        rows.append((product, src, a.GetSymbol(), 1))
        nH = a.GetTotalNumHs()
        if nH:
            rows.append((product, src, "H", nH))

    return (
        pd.DataFrame(rows, columns=["PROD", "EDUCT", "ATOM", "COUNT"])
        .groupby(["PROD", "EDUCT", "ATOM"], as_index=False)["COUNT"]
        .sum()
    )

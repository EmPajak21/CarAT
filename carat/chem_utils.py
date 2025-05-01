from typing import Optional, Union
from rdkit import Chem
import pandas as pd
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.rdchem import Mol
import re


def canonical_smiles(s: Optional[str]) -> Optional[str]:
    """Return the RDKit-canonical SMILES, or None if `s` is None/invalid."""
    if not s:
        return None
    mol = Chem.MolFromSmiles(s)
    return Chem.MolToSmiles(mol) if mol else None


def unmap_smiles(smi: Union[str, Chem.Mol, None]) -> Optional[str]:
    """Strip atom-map numbers from a SMILES string or Mol object."""
    mol = Chem.MolFromSmiles(smi) if isinstance(smi, str) else smi
    if mol:
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        return Chem.MolToSmiles(mol)
    return None


def smiles2molecule(smi: str) -> Union[Mol, None]:
    """Convert smiles to rdkit molecule."""
    return MolFromSmiles(smi) if smi is not None else None


def molecule2smiles(mol: Union[Mol, None]):
    """Convert rdkit molecule to smiles."""
    return MolToSmiles(mol) if mol is not None else None


def contains_carbon(smiles):
    # Use regex to ensure 'C' is matched but not 'Cl'
    if pd.isna(smiles):
        return False
    else:
        return bool(re.search(r"\bC(?!l)\b|C[^l1]", smiles))


def update_smiles(row):
    if (row["SMILES"] == "None") or (row["SMILES"] != row["SMILES"]):
        return "unknown"
    else:
        return canonical_smiles(row["SMILES"])


def bill_of_atoms(mapped_rxn: str) -> pd.DataFrame:
    """Count how many atoms in the *product* originate from each *educt*."""
    left, prod_mapped = mapped_rxn.split(">>")
    educt_frags = left.split(".")

    # Where does each map index come from?
    origin: dict[int, str] = {}
    for frag in educt_frags:
        mol = Chem.MolFromSmiles(frag)
        if not mol:
            continue
        educt = unmap_smiles(frag)
        for a in mol.GetAtoms():
            origin[a.GetAtomMapNum()] = educt

    prod_mol = Chem.MolFromSmiles(prod_mapped)
    product = unmap_smiles(prod_mapped)

    rows: list[tuple[str, str, str, int]] = []
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

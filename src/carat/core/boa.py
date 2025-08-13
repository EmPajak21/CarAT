"""Computes bill-of-atoms (psi) from triplet bills-of-materials.

This module transforms a raw bill-of-materials (BOM) DataFrame into
the psi used in the LP optimisation step. It provides:

- `compute_psi → Dict[str, Dict[Tuple, float]]`: return psi dict
     keyed by triplet, with carbon‐atom shares per reaction leg.
- `make_reaction_smiles`: assemble reaction SMILES strings per triplet.
- `prepare_bom`: filter, normalize and augment a BOM table.

Terminology:
     - BP --> plant or production facility.
     - COCD --> company code.
     - MOV_CAT --> product/reactant.
     - PBG --> material code.
     - PBG_GR --> main product at triplet node.
     - RATIO --> mass ratio of the component.
     - SMILES --> canonical SMILES string of the component.
"""

from typing import Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rxnmapper import RXNMapper

from carat.chem_utils import bill_of_atoms, canonical_smiles


def prepare_bom(bom: pd.DataFrame, *, min_mass: float = 1e-4) -> pd.DataFrame:
    """Subset and normalize the raw bill-of-materials table.

    Parameters
    ----------
    bom : pd.DataFrame
         The raw bill-of-materials DataFrame.

    min_mass : float, optional
         Minimum mass threshold for filtering rows, by default 1e-4.

    Returns
    -------
    pd.DataFrame
         The processed bill-of-materials DataFrame.

    """
    required_columns = {"MOV_CAT", "COCD", "BP", "PBG_GR", "RATIO", "SMILES"}
    missing_columns = required_columns - set(bom.reset_index().columns)
    assert not missing_columns, f"Missing required columns: {missing_columns}"

    bom = bom.copy().reset_index()

    # Mapping MOV_CAT to Product/Reactant
    bom["ROLE"] = bom["MOV_CAT"].map(
        {"GR": "Product", "BY": "Product", "GI": "Reactant"}
    )

    # Creating triplet identifier column
    bom["TRIPLET"] = (
        "t:"
        + bom["COCD"].astype(str)
        + "|"
        + bom["BP"].astype(str)
        + "|"
        + bom["PBG_GR"].astype(str)
    )

    # Randomised AMOUNT placeholder (replace with real data when available)
    bom["AMOUNT"] = np.random.rand(len(bom))
    bom = bom[bom["RATIO"] * bom["AMOUNT"] >= min_mass]

    # Calculating mols to estimate the stoichiometry for reaction SMILES
    bom["AMOUNT"] /= bom.groupby("PBG")["AMOUNT"].transform("sum")
    bom["MASS"] = bom["RATIO"] * bom["AMOUNT"]
    bom["SMILES"] = bom["SMILES"].apply(canonical_smiles)
    bom["MW"] = bom["SMILES"].apply(
        lambda s: ExactMolWt(Chem.MolFromSmiles(s)) if s else np.nan
    )
    bom["MOL"] = bom["MASS"] / bom["MW"]
    return bom


def make_reaction_smiles(bom: pd.DataFrame) -> pd.DataFrame:
    """Create one reaction SMILES per bom produt row.

    Parameters
    ----------
    bom : pd.DataFrame
         The processed bill-of-materials DataFrame.

    Returns
    -------
    pd.DataFrame
         A DataFrame containing reaction SMILES strings for each triplet.

    """
    rows = []
    for trip, grp in bom.groupby("TRIPLET"):
        # Separate reactants and products for the current triplet
        reactants = grp[grp["ROLE"] == "Reactant"]
        products = grp[grp["ROLE"] == "Product"]

        # Iterate over each product in the group
        for _, prod in products.iterrows():
            # Estimate stoichiometry based on ratio of reactant-to-product mols
            n = reactants["MOL"] / prod["MOL"]

            # Rounded and clipped to keep reasonable size of reaction SMILES
            n = n.round().clip(2, 5).astype(int)

            # Construct the left-hand side of the reaction (reactants)
            left = ".".join(
                ".".join([smi] * k)
                for smi, k in zip(reactants["SMILES"], n, strict=False)
            )

            # Append the reaction SMILES (reactants >> product) to the rows list
            rows.append({"TRIPLET": trip, "RXN": f"{left}>>{prod['SMILES']}"})

    # Return a DataFrame containing the reaction SMILES for each triplet
    return pd.DataFrame(rows)


def compute_psi(bom: pd.DataFrame) -> dict[str, dict[tuple[Any, ...], float]]:
    """Return the bill-of-atoms, psi, dict keyed by triplet.

    Parameters
    ----------
    bom : pd.DataFrame
         The raw bill-of-materials DataFrame.

    Returns
    -------
    Dict[str, Dict[Tuple[Any, ...], float]]
         A dictionary keyed by triplet, containing carbon-atom shares per reaction leg.

    """
    bom = prepare_bom(bom)
    rxn_df = make_reaction_smiles(bom)
    mapper = RXNMapper()

    psi = {}
    for triplet, sub in rxn_df.groupby("TRIPLET"):
        boas: list[pd.DataFrame] = []
        for rsmi in sub["RXN"]:
            mapped = mapper.get_attention_guided_atom_maps([rsmi])[0]["mapped_rxn"]
            boas.append(bill_of_atoms(mapped))
        if not boas:
            continue

        boa = pd.concat(boas)

        # Light joins to get PBG + MASS info
        pbg_map = bom.set_index("SMILES")["PBG"].to_dict()
        mass_map = bom.set_index("SMILES")["MASS"].to_dict()
        boa["PROD_PBG"] = boa["PROD"].map(pbg_map)
        boa["EDUCT_PBG"] = boa["EDUCT"].map(pbg_map)
        boa["EDUCT_MASS"] = boa["EDUCT"].map(mass_map)

        # Shares
        grp_cols = ["PROD_PBG", "EDUCT", "ATOM"]
        boa["SMILES_SHARE"] = boa["EDUCT_MASS"] / boa.groupby(grp_cols)[
            "EDUCT_MASS"
        ].transform("sum")
        boa["ATOM_SHARE"] = boa["COUNT"] * boa["SMILES_SHARE"]
        boa["ATOM_SHARE"] /= boa.groupby(["PROD_PBG", "ATOM"])["ATOM_SHARE"].transform(
            "sum"
        )

        psi[triplet] = {
            (
                row["EDUCT_PBG"],
                row["EDUCT"],
                row["PROD_PBG"],
                row["PROD"],
                row["ATOM"],
            ): row["ATOM_SHARE"]
            for _, row in boa[boa["ATOM"] == "C"].iterrows()
        }
    return psi

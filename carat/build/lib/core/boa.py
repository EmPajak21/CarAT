"""Bill-of-atoms ψ calculation from triplet bills‑of‑materials.

Converts bill‑of‑materials DataFrame into the `psi` lookup used downstream
in the linear program optimisation workflow.
"""

from __future__ import annotations

from typing import Dict, Tuple, Any

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rxnmapper import RXNMapper

from chem_utils import canonical_smiles, bill_of_atoms

# ─────────────────────────────────────────────────────────────────────────────
# Pipeline steps
# ─────────────────────────────────────────────────────────────────────────────


def prepare_bom(bom: pd.DataFrame, *, min_mass: float = 1e-4) -> pd.DataFrame:
    """Subset & normalise the raw bill‑of‑materials table."""
    bom = bom.copy().reset_index()
    bom["ROLE"] = bom["MOV_CAT"].map(
        {"GR": "Product", "BY": "Product", "GI": "Reactant"}
    )
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
    bom = bom[bom["RATIO"] * bom["AMOUNT"] >= min_mass]  # drop traces

    bom["AMOUNT"] /= bom.groupby("PBG")["AMOUNT"].transform("sum")
    bom["MASS"] = bom["RATIO"] * bom["AMOUNT"]

    bom["SMILES"] = bom["SMILES"].apply(canonical_smiles)
    bom["MW"] = bom["SMILES"].apply(
        lambda s: ExactMolWt(Chem.MolFromSmiles(s)) if s else np.nan
    )
    bom["MOL"] = bom["MASS"] / bom["MW"]
    return bom


def make_reaction_smiles(bom: pd.DataFrame) -> pd.DataFrame:
    """Create one reaction SMILES per *product* row."""
    rows: list[dict[str, Any]] = []
    for trip, grp in bom.groupby("TRIPLET"):
        reactants = grp[grp["ROLE"] == "Reactant"]
        products = grp[grp["ROLE"] == "Product"]
        for _, prod in products.iterrows():
            # crude stoichiometry: ratio of mols rounded to 1–5
            n = (reactants["MOL"] / prod["MOL"]).round().clip(2, 5).astype(int)
            left = ".".join(
                ".".join([smi] * k) for smi, k in zip(reactants["SMILES"], n)
            )
            rows.append({"TRIPLET": trip, "RXN": f"{left}>>{prod['SMILES']}"})
    return pd.DataFrame(rows)


def compute_psi(bom: pd.DataFrame) -> Dict[str, Dict[Tuple[Any, ...], float]]:
    """Return the ψ‑dictionary keyed by triplet."""
    bom = prepare_bom(bom)
    rxn_df = make_reaction_smiles(bom)
    mapper = RXNMapper()

    psi: Dict[str, Dict[Tuple[Any, ...], float]] = {}
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

        # shares
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

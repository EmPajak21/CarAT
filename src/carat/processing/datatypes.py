"""Data type definitions for CarAT.

Defines custom data types and configurations used in data processing:
- `PostProcessConfig`: Configuration for post-processing tasks.
- `PreProcessConfig`: Configuration for pre-processing tasks.
"""

from dataclasses import dataclass

import networkx as nx
import pandas as pd


@dataclass
class PreProcessConfig:
    """Configuration parameters for initial data preprocessing.

    Attributes
    ----------
    trip_bom : pd.DataFrame
        Triplet bill‐of‐materials DataFrame, indexed by
        (COCD, BP, PBG_GR, PBG, SMILES, ELEMENT).
    dup : pd.DataFrame
        Duplet DataFrame, indexed by (COCD, PBG, SMILES, ELEMENT).
    duplets : Set[Tuple[str, str]]
        Set of all duplet keys (component code, product code).
    inlet_duplets : Set[Tuple[str, str]]
        Subset of duplets considered as inlets.
    var_duplets : Set[Tuple[str, str]]
        Subset of duplets whose attribute mixes are variable.
    triplets : Set[Tuple[str, str, str]]
        Set of all triplet keys (component code, binding point, graph).
    psi : Dict
        Mapping of triplet keys to their ψ share dictionaries.
    nodes : List[str]
        List of graph node identifiers.
    graph : nx.DiGraph
        Directed graph representing the value chain.
    inlets : str, default="base_case_example"
        Type of inlet configuration to apply.

    """

    trip_bom: pd.DataFrame
    dup: pd.DataFrame
    duplets: set[tuple[str, str]]
    inlet_duplets: set[tuple[str, str]]
    var_duplets: set[tuple[str, str]]
    triplets: set[tuple[str, str, str]]
    psi: dict
    nodes: list[str]
    graph: nx.DiGraph
    inlets: str = "base_case_example"


@dataclass
class PostProcessConfig:
    """Processed data ready for LP formulation.

    Attributes
    ----------
    duplets : Set[Tuple[str, str]]
        Set of all duplet keys (component code, product code).
    inlet_duplets : Set[Tuple[str, str]]
        Subset of duplets considered as inlets.
    var_duplets : Set[Tuple[str, str]]
        Subset of duplets whose attribute mixes are variable.
    triplets : Set[Tuple[str, str, str]]
        Set of all triplet keys (component code, binding point, graph).
    trip_bom : pd.DataFrame
        Triplet bill‐of‐materials DataFrame.
    psi : Dict
        Mapping of triplet keys to their psi share dictionaries.
    graph : nx.DiGraph
        Directed graph representing the value chain.
    inlets : str
        Type of inlet configuration applied.
    trip_out : pd.DataFrame
        Filtered triplet output ratios as a DataFrame.
    cps_tank : Set[Tuple[str, str, str]]
        Set of virtual tank configurations (component, product, stream).
    mu : Dict[str, float]
        Normalized connection‐mix share dictionaries for each cps tank.

    """

    duplets: set[tuple[str, str]]
    inlet_duplets: set[tuple[str, str]]
    var_duplets: set[tuple[str, str]]
    triplets: set[tuple[str, str, str]]
    trip_bom: pd.DataFrame
    psi: dict
    graph: nx.DiGraph
    inlets: str

    # New variables from DataPreprocessor
    trip_out: pd.DataFrame
    cps_tank: set[tuple[str, str, str]]
    mu: dict[str, float]

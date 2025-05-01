import pandas as pd
from typing import Dict, List, Set, Tuple, Any
from dataclasses import dataclass


@dataclass
class PreProcessConfig:
    """Configuration parameters for data preprocessing."""

    trip_bom: pd.DataFrame
    dup: pd.DataFrame
    duplets: Set[Tuple[str, str]]
    inlet_duplets: Set[Tuple[str, str]]
    var_duplets: Set[Tuple[str, str]]
    triplets: Set[Tuple[str, str, str]]
    psi: Dict
    nodes: List[str]
    graph: Any
    inlets: str = "base_case_example"


@dataclass
class PostProcessConfig:
    """Processed data output from DataPreprocessor to be used by LPFormulator."""

    # Original fields carried over from DataConfig
    duplets: Set[Tuple[str, str]]
    inlet_duplets: Set[Tuple[str, str]]
    var_duplets: Set[Tuple[str, str]]
    triplets: Set[Tuple[str, str, str]]
    trip_bom: pd.DataFrame
    psi: Dict
    graph: Any
    inlets: str

    # New fields generated during preprocessing
    trip_out: pd.DataFrame
    cps_tank: Any  # Replace with appropriate type
    mu: Dict[str, float]  # Replace with appropriate type
    txt_dict: Dict[str, str]  # Replace with appropriate type

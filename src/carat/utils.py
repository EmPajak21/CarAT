"""General utilities for CarAT.

Helper functions for data parsing, result formatting, and serialization:

- get_example_data(file_path: str = "data/tdi_anon.pkl") -> Dict[str, Any]
- extract_variable_values(var_str: str, model: Any) -> List[Dict[Tuple, float]]
- flatten_dict(data: List[Dict]) -> pd.DataFrame
- format_res_df(var_prefix: str, model: Any) -> pd.DataFrame
- process_results(model: Any) -> Dict[str, pd.DataFrame]
- save_results(results: Dict[str, Any], pkl_output_path: str) -> None
"""

import ast
import pickle
from typing import Any

import pandas as pd


def get_example_data(file_path="data/tdi_anon.pkl"):
    """Load example TDI value chain data from a pickle file and add graph nodes.

    Parameters
    ----------
    file_path : str, optional
        Path to the example data pickle file (default is "data/tdi_anon.pkl").

    Returns
    -------
    Dict[str, Any]
        Dictionary of unpacked data items, with an added "nodes" key containing
        the union of duplets and triplets.

    """
    with open(file_path, "rb") as file:
        data = pickle.load(file)

    # Reformat graph nodes set
    nodes = data["duplets"] | data["triplets"]
    data["nodes"] = nodes
    return data


def extract_variable_values(var_prefix, model):
    """Extract optimization variable values into a list of single-key dicts.

    Parameters
    ----------
    var_prefix : str
        Prefix of the variable names to filter on (e.g., "x_", "z_t").
    model : Any
        Optimized model object with a `.vars` attribute iterable of variable objects
        having `.name` and `.x` attributes.

    Returns
    -------
    List[Dict[Tuple, float]]
        List of dictionaries, each mapping a parsed tuple (derived from var name)
        to its corresponding value.

    """
    res_holder = []
    for v in model.vars:
        if v.name.startswith(var_prefix):
            tuple_str = v.name.split(":", 1)[1].strip()
            res_tuple = ast.literal_eval(tuple_str)
            res_holder.append({res_tuple: v.x})
    return res_holder


def flatten_dict(data):
    """Flatten a list of single-key dictionaries with tuple keys into a DataFrame.

    Parameters
    ----------
    data : List[Dict]
        List where each dict has exactly one key, which may be a tuple.

    Returns
    -------
    pd.DataFrame
        DataFrame where each tuple key is expanded into separate "key_i" columns,
        plus a "value" column for the original dictionary’s value.

    """
    flattened_data = []
    for d in data:
        new_dict = {}
        for k, v in d.items():
            if isinstance(k, tuple):
                # Add each sub-key as its own entry in the new_dict
                for i, sub_key in enumerate(k):
                    new_dict[f"key_{i}"] = sub_key
                # Add the value associated with the tuple key
                new_dict["value"] = v
            else:
                new_dict[k] = v
        flattened_data.append(new_dict)
    return pd.DataFrame(flattened_data)


def format_res_df(var_prefix: str, model: Any) -> pd.DataFrame:
    """Extract decision-variable results from a MIP model into a tidy DataFrame.

    Parameters
    ----------
    var_prefix : str
        Prefix of variable names to include (e.g., "beta_t", "z_d").
    model : Any
        Optimized model object with a `.vars` attribute iterable of variable objects
        having `.name` and `.x` attributes.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns "key_0", "key_1", …, "key_n", and "value",
        where each row corresponds to one decision variable.

    """
    records = []

    for var in getattr(model, "vars", []):
        name_prefix, _, key_part = var.name.partition(":")
        if not name_prefix.startswith(var_prefix):
            continue

        # Safely parse key_part into a tuple by rewrapping it in parentheses
        key_tuple = ast.literal_eval("(" + key_part.strip().strip("()") + ")")

        records.append(list(key_tuple) + [var.x])

    if not records:
        return pd.DataFrame()

    n_keys = len(records[0]) - 1
    columns = [f"key_{i}" for i in range(n_keys)] + ["value"]
    return pd.DataFrame(records, columns=columns)


def process_results(model: Any) -> dict[str, pd.DataFrame]:
    """Process model results into labeled DataFrames for triplet and duplet variables.

    Parameters
    ----------
    model : Any
        Optimized MIP model instance, with decision variables accessible via `.vars`.

    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary containing four DataFrames:
        - "beta_t_res", "beta_d_res": attribute-mix variables for triplets and duplets.
        - "z_t_res", "z_d_res": slack variables for triplets and duplets.

    """
    # Slack variables (positive and negative) for triplets
    z_t_res = format_res_df("z_t", model)
    z_t_res.rename(
        columns={
            "key_0": "COCD",
            "key_1": "BP",
            "key_2": "PBG_GR",
            "key_3": "PBG",
            "key_4": "SMILES",
            "key_5": "ELEMENT",
            "value": "OPTIMAL",
        },
        inplace=True,
    )

    # Slack variables for duplets
    z_d_res = format_res_df("z_d", model)
    z_d_res.rename(
        columns={
            "key_0": "COCD",
            "key_1": "PBG",
            "key_2": "SMILES",
            "key_3": "ELEMENT",
            "value": "OPTIMAL",
        },
        inplace=True,
    )

    # Attribute mix variables for triplets
    beta_t_res = format_res_df("beta_t", model)
    beta_t_res.rename(
        columns={
            "key_0": "COCD",
            "key_1": "BP",
            "key_2": "PBG_GR",
            "key_3": "PBG",
            "key_4": "SMILES",
            "key_5": "ELEMENT",
            "key_6": "ATTRIBUTE",
            "value": "OPTIMAL",
        },
        inplace=True,
    )

    # Attribute mix variables for duplets
    beta_d_res = format_res_df("beta_d", model)
    beta_d_res.rename(
        columns={
            "key_0": "COCD",
            "key_1": "PBG",
            "key_2": "SMILES",
            "key_3": "ELEMENT",
            "key_4": "ATTRIBUTE",
            "value": "OPTIMAL",
        },
        inplace=True,
    )

    return {
        "beta_t_res": beta_t_res,
        "beta_d_res": beta_d_res,
        "z_t_res": z_t_res,
        "z_d_res": z_d_res,
    }


def save_results(results: dict[str, Any], pkl_output_path: str) -> None:
    """Save a results dictionary to a pickle file.

    Parameters
    ----------
    results : Dict[str, Any]
        Dictionary containing result DataFrames and other serializable objects.
    pkl_output_path : str
        File path where the pickle file will be written.

    """
    with open(pkl_output_path, "wb") as f:
        pickle.dump(results, f)

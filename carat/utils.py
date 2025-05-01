import pickle
import ast
import pandas as pd
from typing import Dict, Any
import ast


def get_example_data(file_path="data/tdi_anon.pkl"):
    """
    Load example TDI value chain data from a pickle file.

    Parameters:
    file_path : str
        Path to the pickle file containing the example data.

    Returns:
    tuple : Unpacked data items and nodes
    """

    with open(file_path, "rb") as file:
        data = pickle.load(file)

    # Reformat graph nodes set
    nodes = data["duplets"] | data["triplets"]
    data["nodes"] = nodes
    return data


def opt_res_to_list_dict(var_str, model):
    res_holder = []
    for v in model.vars:
        if v.name.startswith(var_str):
            tuple_str = v.name.split(":", 1)[1].strip()
            res_tuple = ast.literal_eval(tuple_str)
            res_holder.append({res_tuple: v.x})
    return res_holder


def flatten_dict(data):
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
    """
    Extract decision-variable results from a MIP model into a tidy DataFrame.
    """
    records = []

    for var in getattr(model, "vars", []):
        name_prefix, _, key_part = var.name.partition(":")
        if not name_prefix.startswith(var_prefix):
            continue

        # strip surrounding parentheses then re-add them so ast can parse
        key_tuple = ast.literal_eval("(" + key_part.strip().strip("()") + ")")
        # key_tuple is now a real Python tuple like (0, 1, 2)

        records.append(list(key_tuple) + [var.x])

    if not records:
        return pd.DataFrame()

    n_keys = len(records[0]) - 1
    columns = [f"key_{i}" for i in range(n_keys)] + ["value"]
    return pd.DataFrame(records, columns=columns)


def process_results(model: Any, txt_dict: Dict[str, Any]) -> Dict[str, pd.DataFrame]:
    """
    Process model results into labeled DataFrames for triplet and duplet variables.

    Args:
        model: Optimized MIP model instance.
        txt_dict: Dictionary for text representations (returned but not modified).

    Returns:
        Dictionary with DataFrames: beta_t_res, beta_d_res, z_t_res, z_d_res.
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


def save_results(results: Dict[str, Any], pkl_output_path: str) -> None:
    """
    Save results dictionary to a pickle file.

    Args:
        results: Dictionary containing result DataFrames and other objects.
        pkl_output_path: Path to write the pickle file.
    """
    with open(pkl_output_path, "wb") as f:
        pickle.dump(results, f)

"""Formulates and solves an LP for biogenic carbon content calculation.

This module takes the output of a DataPreprocessor (a preprocessed value
chain) and builds a Linear Program model (via Python-MIP/CBC) to trace
and optimize the share of biogenic vs. fossil carbon throughout the network.
Results can be processed and serialized via carat.utils.
"""

import logging
from typing import Any

from mip import CBC, Model, minimize, xsum

from carat.utils import process_results, save_results

# Configure logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class LPFormulator:
    """Formulate and solve linear programming optimization for value chain analysis.

    Parameters
    ----------
    preprocessed_data : Dict[str, Any]
        Preprocessed data from the DataPreprocessor.
    inlets : str, optional
        Inlet conditions to use, by default "base_case_example".

    """

    def __init__(
        self, preprocessed_data: dict[str, Any], inlets: str = "base_case_example"
    ):
        """Initialize the LP formulator with preprocessed data.

        Args:
            preprocessed_data: Dictionary containing preprocessed value chain data
            inlets: Type of inlet configuration to use (default: "base_case_example")

        """
        for key, value in vars(preprocessed_data).items():
            setattr(self, key, value)

        self.model = Model(solver_name=CBC)
        self.A, self.e = self._define_sets()
        self.decision_vars = self._define_decision_vars()
        self.inlets = inlets

    def solve(self, output_results: bool = True, pkl_output_path: str = "") -> dict:
        """Construct and solve the LP optimization problem.

        Parameters
        ----------
        output_results : bool, optional
            Whether to output results to the logger, by default True.
        pkl_output_path : str, optional
            Path to save results as a pickle file, by default "".

        Returns
        -------
        Dict
            Results of the optimization.

        """
        self._set_inlet_conditions()
        self._add_attribute_sum_constraints()
        self._add_triplet_beta_constraints()
        self._add_duplet_beta_constraints()
        self._set_objective_function()

        self.model.solver_options = [("tol", 1e-4)]
        self.model.optimize()

        if output_results:
            results = process_results(self.model)
            logger.info("Beta_d:\n%s", results["beta_d_res"])
            logger.info("Beta_t:\n%s", results["beta_t_res"])
            logger.info("z_d:\n%s", results["z_d_res"])
            logger.info("z_t:\n%s", results["z_t_res"])

        if pkl_output_path:
            save_results(
                results,
                pkl_output_path,
            )

        return results

    def _define_sets(self) -> tuple[set[str], str]:
        """Define key sets for the LP formulation.

        Returns
        -------
        Tuple[Set[str], str]
            A set of attributes and the element to be traced.

        """
        A = {"biogenic", "fossil"}
        e = "C"
        return A, e

    def _define_decision_vars(self) -> dict[str, dict]:
        """Define all decision variables for the LP model.

        Returns
        -------
        Dict[str, Dict]
            Dictionary containing decision variables for duplets and triplets.

        """
        beta_d = {}
        z_d_pos = {}
        z_d_neg = {}

        for c, p, s in self.cps_tank:
            z_d_pos[c, p, s, self.e] = self.model.add_var(
                name=f"z_d_pos:{c, p, s, self.e}", lb=0, ub=1
            )
            z_d_neg[c, p, s, self.e] = self.model.add_var(
                name=f"z_d_neg:{c, p, s, self.e}", lb=-1, ub=0
            )
            for a in self.A:
                beta_d[c, p, s, self.e, a] = self.model.add_var(
                    name=f"beta_d:{c, p, s, self.e, a}", lb=0, ub=1
                )

        beta_t = {}
        z_t_pos = {}
        z_t_neg = {}
        for c, b, g, p, s in self.trip_out:
            z_t_pos[c, b, g, p, s, self.e] = self.model.add_var(
                name=f"z_t_pos:{c, b, g, p, s, self.e}", lb=0, ub=1
            )
            z_t_neg[c, b, g, p, s, self.e] = self.model.add_var(
                name=f"z_t_neg:{c, b, g, p, s, self.e}", lb=-1, ub=0
            )
            for a in self.A:
                beta_t[c, b, g, p, s, self.e, a] = self.model.add_var(
                    name=f"beta_t:{c, b, g, p, s, self.e, a}", lb=0, ub=1
                )

        return {
            "beta_d": beta_d,
            "z_d_pos": z_d_pos,
            "z_d_neg": z_d_neg,
            "beta_t": beta_t,
            "z_t_pos": z_t_pos,
            "z_t_neg": z_t_neg,
        }

    def _set_inlet_conditions(self):
        """Set inlet conditions based on the selected inlet type.

        - If `self.inlets` is `"base_case_example"`, all inlet duplets are forced to
        fossil = 1 and biogenic = 0.
        - If `self.inlets` is `"C1"`, the specified inlet node (`("COMP_2", "PROD_12")`)
        is switched to biogenic = 1 (and fossil = 0), while all other inlet duplets
        remain fossil = 1 and biogenic = 0.
        """
        beta_d = self.decision_vars["beta_d"]

        if self.inlets == "base_case_example":
            for c, p, s, e, _a in beta_d:
                if (c, p) in self.inlet_duplets:
                    self.model += beta_d[c, p, s, e, "fossil"] == 1
                    self.model += beta_d[c, p, s, e, "biogenic"] == 0

        elif self.inlets == "C1":
            inlet_node = ("COMP_2", "PROD_12")
            for c, p, s, e, _a in beta_d:
                if (c, p) in self.inlet_duplets:
                    if c == inlet_node[0] and p == inlet_node[1]:
                        logger.debug("Switching to biogenic for inlet %s", inlet_node)
                        self.model += beta_d[c, p, s, e, "fossil"] == 0
                        self.model += beta_d[c, p, s, e, "biogenic"] == 1
                    else:
                        self.model += beta_d[c, p, s, e, "fossil"] == 1
                        self.model += beta_d[c, p, s, e, "biogenic"] == 0

    def _add_attribute_sum_constraints(self):
        """Add constraints ensuring attribute shares sum to 1.

        For each triplet (β_t) and each duplet (β_d), enforce that the
        sum of the attribute‐share variables minus their slack variables equals 1.
        """
        beta_d = self.decision_vars["beta_d"]
        z_d_pos = self.decision_vars["z_d_pos"]
        z_d_neg = self.decision_vars["z_d_neg"]
        beta_t = self.decision_vars["beta_t"]
        z_t_pos = self.decision_vars["z_t_pos"]
        z_t_neg = self.decision_vars["z_t_neg"]

        for c, b, g, p, s, e, _ in beta_t:
            self.model += (
                xsum(beta_t[c, b, g, p, s, e, a] for a in self.A)
                - z_t_pos[c, b, g, p, s, e]
                - z_t_neg[c, b, g, p, s, e]
                == 1
            )

        for c, p, s, e, _ in beta_d:
            self.model += (
                xsum(beta_d[c, p, s, e, a] for a in self.A)
                - z_d_pos[c, p, s, e]
                - z_d_neg[c, p, s, e]
                == 1
            )

    def _add_triplet_beta_constraints(self):
        """Add elemental attribute constraints for triplets.

        This method builds constraints equating each β_t variable to the
        weighted sum of β_d variables via the psi coefficients.
        """
        beta_d = self.decision_vars["beta_d"]
        beta_t = self.decision_vars["beta_t"]

        for c, b, g in self.triplets:
            trip_key = f"t:{c}|{b}|{g}"
            psi_t = self.psi[trip_key]

            beta_t_trip = {
                key: var
                for key, var in beta_t.items()
                if key[0] == c and key[1] == b and key[2] == g
            }

            for (c_t, _b_t, _g_t, p_t, s_t, _e_t, a_t), var in beta_t_trip.items():
                cumulative = 0
                for (p, s, p_pr, s_pr, e), val in psi_t.items():
                    if p_pr == p_t and s_pr == s_t:
                        cumulative += val * beta_d[c_t, p, s, e, a_t]
                self.model += var == cumulative

    def _add_duplet_beta_constraints(self):
        """Add elemental attribute constraints for duplets.

        For each duplet in `cps_tank` that’s a variable, this method:
        - Gathers matching triplet β_t variables for the same component,
            binding point, graph, stream, and attribute.
        - Renormalizes their ψ shares if they don’t sum to 1.
        - Builds a weighted sum of those β_t variables.
        - Constrains each duplet β_d to equal that weighted sum.
        """
        beta_d = self.decision_vars["beta_d"]
        beta_t = self.decision_vars["beta_t"]

        for c, p, s in self.cps_tank:
            if (c, p) in self.var_duplets:
                mu_cps = self.mu[(c, p, s)]
                for a in self.A:
                    cumulative = 0
                    shares = []
                    terms = []
                    for (c_t, b_t, g_t), share in mu_cps.items():
                        matches = [
                            (key, var)
                            for key, var in beta_t.items()
                            if key[0] == c_t
                            and key[1] == b_t
                            and key[2] == g_t
                            and key[4] == s
                            and key[6] == a
                        ]
                        if matches:
                            key_t, var_t = matches[0]
                            shares.append(share)
                            terms.append(var_t)
                    if shares and sum(shares) != 1:
                        norm = [sh / sum(shares) for sh in shares]
                        cumulative = sum(
                            n * t for n, t in zip(norm, terms, strict=False)
                        )
                    else:
                        cumulative = sum(
                            share * terms[i] for i, share in enumerate(shares)
                        )
                    if shares:
                        self.model += beta_d[c, p, s, self.e, a] == cumulative

    def _set_objective_function(self):
        """Set the objective function to minimize total slack.

        Model objective is the absolute sum of slack variables
        penalizing deviations in both duplet and triplet constraints.
        """
        z_d_pos = self.decision_vars["z_d_pos"]
        z_d_neg = self.decision_vars["z_d_neg"]
        z_t_pos = self.decision_vars["z_t_pos"]
        z_t_neg = self.decision_vars["z_t_neg"]

        self.model.objective = minimize(
            -xsum(z_d_neg.values())
            - xsum(z_t_neg.values())
            + xsum(z_d_pos.values())
            + xsum(z_t_pos.values())
        )

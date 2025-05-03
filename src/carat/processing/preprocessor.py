"""Data preprocessing utilities for CarAT.

Main class:
- DataPreprocessor: transform raw value-chain configurations into
structured data for LP formulation.
"""

from typing import Dict, Set, Tuple

from carat.chem_utils import contains_carbon
from carat.processing.datatypes import PostProcessConfig, PreProcessConfig


class DataPreprocessor:
    """
    Preprocess data for linear programming formulation of value chain analysis.

    Parameters
    ----------
    config : PreProcessConfig
        Configuration parameters for initial data preprocessing.
    """

    def __init__(self, config: PreProcessConfig):
        """
        Initialize the preprocessor and prepare all derived data.

        Parameters
        ----------
        config : PreProcessConfig
            Input configuration holding raw data and settings.
        """
        # Unpack config into individual instance attributes
        for key, value in vars(config).items():
            setattr(self, key, value)

        # Initialize additional attributes
        self._prepare_derived_data()

    def preprocess(self) -> PostProcessConfig:
        """
        Return processed data ready for LP formulation.

        Returns
        -------
        PostProcessConfig
            Dataclass containing both original and derived datasets.
        """
        # Preprocessing is handled in the constructor; return processed data.
        return PostProcessConfig(
            duplets=self.duplets,
            inlet_duplets=self.inlet_duplets,
            var_duplets=self.var_duplets,
            triplets=self.triplets,
            trip_bom=self.trip_bom,
            psi=self.psi,
            graph=self.graph,
            inlets=self.inlets,
            trip_out=self.trip_out,
            cps_tank=self.cps_tank,
            mu=self.mu,
        )

    def _prepare_derived_data(self):
        """
        Prepare all derived data required for LP formulation.

        This method computes:
        - trip_out: Filtered triplet outputs.
        - cps_tank: Virtual tank configurations.
        - mu: Connection mix shares.
        """
        self.trip_out = self._get_trip_out()
        self.cps_tank = self._get_cps_tank()
        self.mu = self._get_con_mix_share()

    def _get_trip_out(self) -> Dict:
        """
        Extract and process triplet output ratios.

        Filters for product entries, carbon‐containing streams,
        and reorders index levels.

        Returns
        -------
        Dict[Tuple[str, str, str, str, str], float]
            Mapping from triplet keys to their output ratio.
        """

        # This returns a DataFrame of all product entries with carbon atoms.
        trip_out = self.trip_bom.loc[
            # Filter out rows where MOV_CAT is "GI" i.e., reactants and SMILES.
            (self.trip_bom["MOV_CAT"] != "GI")
            & (
                # Includ SMILES with carbon atoms.
                self.trip_bom.index.get_level_values(4).str.contains("C")
                | self.trip_bom.index.get_level_values(4).str.contains("unknown")
            ),
            "RATIO",  # Select the "RATIO" column
        ]

        col_names = ["COCD", "BP", "PBG_GR", "PBG", "SMILES"]
        trip_out = trip_out.reorder_levels(col_names)
        trip_out = trip_out[
            trip_out.index.get_level_values("SMILES").map(contains_carbon)
        ]

        # Assert that the names are exactly as expected

        trip_out = trip_out.reorder_levels(col_names)
        assert (
            trip_out.index.names == col_names
        ), f"trip_out.index.names = {trip_out.index.names!r}, expected {col_names!r}"
        return trip_out.to_dict()

    def _get_cps_tank(self) -> Set[Tuple[str, str, str]]:
        """
        Extract virtual tank configurations for duplets.

        Returns
        -------
        Set[Tuple[str, str, str]]
            Each tuple is (component, product, stream) for tanks containing carbon.
        """
        col_names = ["COCD", "PBG", "SMILES"]
        cps_tank = self.dup.reset_index()[col_names].drop_duplicates()
        cps_tank = cps_tank[cps_tank["SMILES"].apply(contains_carbon)]
        return {tuple(row) for row in cps_tank.itertuples(index=False)}

    def _get_con_mix_share(self, renorm_mu: bool = True) -> Dict:
        """
        Calculate connection‐mix share for each cps tank.

        Parameters
        ----------
        renorm_mu : bool, optional
            If True, renormalize shares so they sum to 1 when needed (default is True).

        Returns
        -------
        Dict[Tuple[str, str, str], Dict[Tuple, float]]
            Maps each (component, product, stream) to a dict of share values.
        """
        mu = {}
        for c, p, s in self.cps_tank:
            dup_cp = self.dup.xs((c, p, s), level=("COCD", "PBG", "SMILES"))[
                "SHARE"
            ].to_dict()
            if (c, p) in self.var_duplets and renorm_mu:
                total = sum(dup_cp.values())
                if total != 1:
                    dup_cp = {k: v / total for k, v in dup_cp.items()}
            mu[(c, p, s)] = dup_cp
        return mu

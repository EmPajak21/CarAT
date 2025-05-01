import logging
from typing import Dict, Set, Tuple
from carat.chem_utils import contains_carbon

from carat.processing.datatypes import PreProcessConfig, PostProcessConfig

# Configure logger
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class DataPreprocessor:
    """Preprocess data for linear programming formulation of value chain analysis."""

    def __init__(self, config: PreProcessConfig):
        # Unpack config into individual instance attributes
        for key, value in vars(config).items():
            setattr(self, key, value)

        # Initialize additional attributes
        self._prepare_derived_data()

    def preprocess(self) -> PostProcessConfig:
        """Main method to preprocess data and return processed data."""
        # The preprocessing logic is already encapsulated in the constructor
        # so we can just return the processed data here.
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
            txt_dict=self.txt_dict,
        )

    def _prepare_derived_data(self):
        """Prepare all derived data required for LP formulation."""
        self.trip_out = self._get_trip_out()
        self.cps_tank = self._get_cps_tank()
        self.mu = self._get_con_mix_share()
        self.txt_dict = self._get_txt_dictionaries()

    def _get_trip_out(self) -> Dict:
        """Extract and process triplet output data."""
        trip_out = self.trip_bom.loc[
            (self.trip_bom["MOV_CAT"] != "GI")
            & (
                self.trip_bom.index.get_level_values(4).str.contains("C")
                | self.trip_bom.index.get_level_values(4).str.contains("unknown")
            ),
            "RATIO",
        ]
        trip_out = trip_out.reorder_levels(["COCD", "BP", "PBG_GR", "PBG", "SMILES"])
        trip_out = trip_out[
            trip_out.index.get_level_values("SMILES").map(contains_carbon)
        ]
        return trip_out.to_dict()

    def _get_cps_tank(self) -> Set[Tuple[str, str, str]]:
        """Extract virtual tank configurations."""
        cps_tank = self.dup.reset_index()[["COCD", "PBG", "SMILES"]].drop_duplicates()
        cps_tank = cps_tank[cps_tank["SMILES"].apply(contains_carbon)]
        return {tuple(row) for row in cps_tank.itertuples(index=False)}

    def _get_con_mix_share(self, renorm_mu: bool = True) -> Dict:
        """Calculate connection mix share for sub value chain."""
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

    def _get_txt_dictionaries(self) -> Dict:
        """Create dictionaries for text representations."""
        trip_bom = self.trip_bom.reset_index()
        trip_bom["TRIPLET"] = (
            "t:" + trip_bom["COCD"] + "|" + trip_bom["BP"] + "|" + trip_bom["PBG_GR"]
        )
        txt_dict = trip_bom[["TRIPLET", "BP_TXT"]].set_index("TRIPLET").to_dict()

        dup = self.dup.reset_index()
        dup["DUPLET"] = "d:" + dup["COCD"] + "|" + dup["PBG"]
        pbg_txt_dict = dup[["DUPLET", "PBG_TXT"]].set_index("DUPLET").to_dict()
        txt_dict.update(pbg_txt_dict)

        return txt_dict

"""Conservation-audited zero-dimensional plasma-reactor building blocks.

The package separates universal bookkeeping (species, reactions, geometry,
units) from chemistry decks and equipment-specific closures.  A global model
is only as predictive as its rate and wall-loss evidence; these classes make
that evidence and every open-system exchange explicit.
"""

from .geometry import CylindricalReactor, ElectropositiveEdgeFactors
from .argon import (
    ARGON_4S_METASTABLE_ENERGY_EV,
    ARGON_IONIZATION_ENERGY_EV,
    ARGON_METASTABLE_IONIZATION_ENERGY_EV,
    build_lee_lieberman_argon_volume_network,
    lee_lieberman_argon_species,
)
from .network import (
    CM3_TO_M3,
    E_CHARGE_C,
    ConstantRateCoefficient,
    ElectronArrheniusRateCoefficient,
    RateContext,
    Reaction,
    ReactionNetwork,
    Species,
)
from .model import (
    ArgonGlobalCondition,
    ArgonGlobalSolution,
    ArgonTransportProvider,
    ArgonTransportState,
    FixedArgonTransportProvider,
    LeeLiebermanArgonGlobalModel,
    PASCAL_PER_MTORR,
)

__all__ = [
    "CM3_TO_M3",
    "E_CHARGE_C",
    "ARGON_4S_METASTABLE_ENERGY_EV",
    "ARGON_IONIZATION_ENERGY_EV",
    "ARGON_METASTABLE_IONIZATION_ENERGY_EV",
    "ArgonGlobalCondition",
    "ArgonGlobalSolution",
    "ArgonTransportProvider",
    "ArgonTransportState",
    "ConstantRateCoefficient",
    "CylindricalReactor",
    "ElectronArrheniusRateCoefficient",
    "ElectropositiveEdgeFactors",
    "FixedArgonTransportProvider",
    "LeeLiebermanArgonGlobalModel",
    "PASCAL_PER_MTORR",
    "RateContext",
    "Reaction",
    "ReactionNetwork",
    "Species",
    "build_lee_lieberman_argon_volume_network",
    "lee_lieberman_argon_species",
]

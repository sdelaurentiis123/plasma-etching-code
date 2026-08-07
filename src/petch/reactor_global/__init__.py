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
from .transport import (
    ARGON_MASS_KG,
    ATOMIC_MASS_UNIT_KG,
    STANDARD_PRESSURE_PA,
    LeeLiebermanArgonTransportProvider,
    argon_relative_temperature_eV,
    lee_lieberman_argon_ion_temperature_eV,
    nist_argon_self_diffusion_m2_s,
    phelps_argon_momentum_transfer_cross_section_m2,
    phelps_argon_momentum_transfer_rate_m3_s,
)

__all__ = [
    "CM3_TO_M3",
    "E_CHARGE_C",
    "ARGON_4S_METASTABLE_ENERGY_EV",
    "ARGON_IONIZATION_ENERGY_EV",
    "ARGON_METASTABLE_IONIZATION_ENERGY_EV",
    "ARGON_MASS_KG",
    "ATOMIC_MASS_UNIT_KG",
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
    "LeeLiebermanArgonTransportProvider",
    "PASCAL_PER_MTORR",
    "RateContext",
    "Reaction",
    "ReactionNetwork",
    "Species",
    "STANDARD_PRESSURE_PA",
    "argon_relative_temperature_eV",
    "build_lee_lieberman_argon_volume_network",
    "lee_lieberman_argon_species",
    "lee_lieberman_argon_ion_temperature_eV",
    "nist_argon_self_diffusion_m2_s",
    "phelps_argon_momentum_transfer_cross_section_m2",
    "phelps_argon_momentum_transfer_rate_m3_s",
]

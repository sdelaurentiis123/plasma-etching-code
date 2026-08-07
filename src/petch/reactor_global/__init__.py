"""Conservation-audited zero-dimensional plasma-reactor building blocks.

The package separates universal bookkeeping (species, reactions, geometry,
units) from chemistry decks and equipment-specific closures.  A global model
is only as predictive as its rate and wall-loss evidence; these classes make
that evidence and every open-system exchange explicit.
"""

from .geometry import CylindricalReactor, ElectropositiveEdgeFactors
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

__all__ = [
    "CM3_TO_M3",
    "E_CHARGE_C",
    "ConstantRateCoefficient",
    "CylindricalReactor",
    "ElectronArrheniusRateCoefficient",
    "ElectropositiveEdgeFactors",
    "RateContext",
    "Reaction",
    "ReactionNetwork",
    "Species",
]

#!/usr/bin/env python3
"""Compare petch's two-term EEPF with a local LGPL BOLOS installation.

The oracle is optional and never imported by the petch package.  This script
uses only a manufactured, redistributable collision deck and prints a JSON
receipt to stdout.  No BOLSIG+, BOLOS, or LXCat source/data bytes are copied
into the repository.
"""
from __future__ import annotations

import argparse
from importlib import metadata
import json
from pathlib import Path
import sys

import numpy as np
import scipy
import scipy.integrate

from petch.reactor_global.electron_collision_deck import (
    parse_bolsig_lxcat_bytes,
)
from petch.reactor_global.electron_kinetics import (
    DeterministicTwoTermBoltzmannSolver,
    ElectronEnergyGrid,
    TwoTermBoltzmannCondition,
    normalize_eepf,
)


MANUFACTURED_DECK = b"""ELASTIC
manufactured
1.0e-5 / electron-to-target mass ratio
------------------------------------------------------------
0.0 2.0e-20
300.0 2.0e-20
------------------------------------------------------------

EXCITATION
manufactured -> manufactured(v=1)
2.0 / threshold energy
------------------------------------------------------------
0.0 0.0
2.0 0.0
5.0 1.2e-20
300.0 1.2e-20
------------------------------------------------------------
"""


def _relative_error(value: float, reference: float) -> float:
    return abs(value - reference) / abs(reference)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--bolos-path",
        type=Path,
        required=True,
        help="directory containing the locally installed/importable bolos package",
    )
    parser.add_argument(
        "--cells",
        type=int,
        nargs="+",
        default=(600, 1200, 2400),
    )
    args = parser.parse_args()
    sys.path.insert(0, str(args.bolos_path.resolve()))

    # BOLOS 0.2 predates SciPy's simps -> simpson rename.  This compatibility
    # alias changes no quadrature definition and remains local to the audit.
    if not hasattr(scipy.integrate, "simps"):
        scipy.integrate.simps = (
            lambda y, x=None, **kwargs: scipy.integrate.simpson(
                y, x=x, **kwargs))

    from bolos import grid as bolos_grid  # noqa: PLC0415
    from bolos import solver as bolos_solver  # noqa: PLC0415

    deck = parse_bolsig_lxcat_bytes(
        MANUFACTURED_DECK,
        source_database="manufactured-redistributable-oracle-deck",
        retrieved_at="2026-08-08",
        source_reference=__file__,
        target="manufactured",
    )
    maximum_energy_eV = 200.0
    reduced_field_Td = 5.0
    gas_temperature_K = 300.0
    rows = []
    for cell_count in args.cells:
        grid = ElectronEnergyGrid.linear(maximum_energy_eV, cell_count)
        condition = TwoTermBoltzmannCondition(
            reduced_electric_field_Td=reduced_field_Td,
            gas_temperature_K=gas_temperature_K,
            target_mole_fractions={"manufactured": 1.0},
        )
        petch_solution = DeterministicTwoTermBoltzmannSolver(
            grid, deck).solve(
                condition,
                damping=1.0,
                relative_tolerance=1.0e-10,
                maximum_tail_population_fraction=1.0e-6,
            )

        oracle = bolos_solver.BoltzmannSolver(
            bolos_grid.LinearGrid(0.0, maximum_energy_eV, cell_count))
        oracle.add_process(
            kind="ELASTIC",
            target="manufactured",
            mass_ratio=1.0e-5,
            data=np.array([[0.0, 2.0e-20], [300.0, 2.0e-20]]),
        )
        oracle.add_process(
            kind="EXCITATION",
            target="manufactured",
            product="manufactured(v=1)",
            threshold=2.0,
            data=np.array([
                [0.0, 0.0],
                [2.0, 0.0],
                [5.0, 1.2e-20],
                [300.0, 1.2e-20],
            ]),
        )
        oracle.target["manufactured"].density = 1.0
        oracle.kT = (
            1.380649e-23 * gas_temperature_K / 1.602176634e-19)
        oracle.EN = reduced_field_Td * 1.0e-21
        oracle.init()
        oracle_eepf = oracle.converge(
            oracle.maxwell(2.0), rtol=1.0e-10, maxn=300)
        oracle_mean = float(oracle.mean_energy(oracle_eepf))
        oracle_rate = float(oracle.rate(
            oracle_eepf,
            oracle.search("manufactured", "manufactured(v=1)"),
        ))
        oracle_minimum = float(np.min(oracle_eepf))
        oracle_scale = max(float(np.max(oracle_eepf)), 1.0)
        if oracle_minimum < -1.0e-12 * oracle_scale:
            raise RuntimeError("BOLOS oracle returned a physical-scale negative EEPF")
        oracle_negative_population = float(np.dot(
            np.maximum(-oracle_eepf, 0.0), grid.normalization_weights))
        normalized_oracle_eepf = normalize_eepf(
            grid, np.maximum(oracle_eepf, 0.0))
        petch_mean = float(petch_solution.distribution.mean_energy_eV)
        petch_rate = float(
            petch_solution.collision_moments[1].rate_coefficient_m3_s)
        rows.append({
            "cell_count": cell_count,
            "petch_mean_energy_eV": petch_mean,
            "oracle_mean_energy_eV": oracle_mean,
            "mean_energy_relative_error": _relative_error(
                petch_mean, oracle_mean),
            "petch_excitation_rate_m3_s": petch_rate,
            "oracle_excitation_rate_m3_s": oracle_rate,
            "excitation_rate_relative_error": _relative_error(
                petch_rate, oracle_rate),
            "eepf_weighted_l1": float(np.sum(
                np.abs(
                    petch_solution.distribution.eepf_eV_minus_3_over_2
                    - normalized_oracle_eepf
                ) * grid.normalization_weights
            )),
            "oracle_minimum_raw_eepf": oracle_minimum,
            "oracle_roundoff_negative_population_fraction": (
                oracle_negative_population),
        })

    try:
        bolos_version = metadata.version("bolos")
    except metadata.PackageNotFoundError:
        bolos_version = "unknown-local-install"
    receipt = {
        "schema": "petch.two_term_bolos_oracle.v1",
        "manufactured_deck_sha256": deck.payload_sha256,
        "condition": {
            "reduced_electric_field_Td": reduced_field_Td,
            "gas_temperature_K": gas_temperature_K,
            "maximum_energy_eV": maximum_energy_eV,
        },
        "oracle": {
            "name": "BOLOS",
            "version": bolos_version,
            "license": "LGPL; local audit dependency only",
        },
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "rows": rows,
        "finest_grid_pass": (
            rows[-1]["mean_energy_relative_error"] < 0.01
            and rows[-1]["excitation_rate_relative_error"] < 0.01
        ),
        "supports_direct_swarm_grade": False,
        "supports_reactor_state_prediction": False,
        "supports_wafer_flux": False,
        "supports_feature_depth": False,
    }
    print(json.dumps(receipt, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

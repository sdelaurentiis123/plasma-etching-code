#!/usr/bin/env python3
"""Replay Metzler's cyclic C4F8/Ar+ Si and SiO2 surface-state board.

The board deliberately keeps three different observables separate:

* ellipsometric deposited-film thickness and etched depth per cycle;
* XPS-derived F/C within the fluorocarbon film;
* ``delta F/C``, Metzler's difference between two XPS ratios, used by the
  source as a *relative* proxy for fluorine transferred into the substrate.

Neither film thickness nor delta-F/C is silently converted into an absolute
areal atom inventory.  The source does not publish the density/composition
needed for the first conversion or an XPS sensitivity/depth model for the
second.

The official thesis PDF is not redistributed.  This script verifies its
checksum, renders the five audited pages at 240 dpi, checks their pixels,
reproduces the committed CSVs and manifest, and optionally draws PIL QA
overlays at full raster resolution.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tempfile

from PIL import Image, ImageDraw


ROOT = Path(__file__).resolve().parents[1]
OUTPUT_DIR = ROOT / "data" / "experimental" / "metzler_2016"
DEPTH_CSV = OUTPUT_DIR / "figures6_5_6_6_cyclic_depth.csv"
YIELD_CSV = OUTPUT_DIR / "figure6_9_cycle_averaged_yield.csv"
XPS_CSV = OUTPUT_DIR / "figures6_14_6_15_xps_cycle_state.csv"
MANIFEST_PATH = OUTPUT_DIR / "digitization_manifest.json"
DEFAULT_PDF = ROOT / "tmp" / "sources" / "metzler_2016" / "thesis.pdf"

SOURCE = {
    "citation": (
        "D. Metzler, High Precision Plasma Etch for Pattern Transfer: "
        "Towards Fluorocarbon Based Atomic Layer Etching, PhD thesis, "
        "University of Maryland (2016)"
    ),
    "url": (
        "https://api.drum.lib.umd.edu/server/api/core/bitstreams/"
        "227fdb28-6ea7-4a8d-92ae-5c521f2d1e0b/content"
    ),
    "item_url": (
        "https://drum.lib.umd.edu/items/"
        "b01bb413-6152-492a-98ec-0abe752ea240"
    ),
    "pdf_sha256": (
        "ea5701d0bcf67b56403625253f7bb619c0e3e3a5e0a9cfd2ff6f1e435fa90f62"
    ),
    "pages": {
        141: {
            "figure": "6.5",
            "png_sha256": (
                "65533648fdf3047aa42c15cc4fffabe394eb6577f7320cf70a187dd2e4f5b898"
            ),
        },
        142: {
            "figure": "6.6",
            "png_sha256": (
                "52845fe4c47ca9ba8f270af1a0de4fc71a04cc489849a51044cf1c5cc4ebc5fa"
            ),
        },
        147: {
            "figure": "6.9",
            "png_sha256": (
                "86e665a6d1fc417ac9ae6ff8b7b900e36f69871a5b84d0467756b7cf2d5824e2"
            ),
        },
        157: {
            "figure": "6.14",
            "png_sha256": (
                "10ecc9e1693461ae6c75d98b554ccf42e5bacb0194f3dfab824a94ef96393e33"
            ),
        },
        158: {
            "figure": "6.15",
            "png_sha256": (
                "db6707764a20786f033c5f4ff40d90234502d663b1fe637bcc8415f41957437e"
            ),
        },
    },
}


@dataclass(frozen=True)
class Axis:
    x_left_px: float
    x_right_px: float
    x_left_value: float
    x_right_value: float
    y_zero_px: float
    y_reference_px: float
    y_zero_value: float
    y_reference_value: float

    def x_value(self, x_px: float) -> float:
        fraction = (
            (float(x_px) - self.x_left_px)
            / (self.x_right_px - self.x_left_px)
        )
        return self.x_left_value + fraction * (
            self.x_right_value - self.x_left_value
        )

    def y_value(self, y_px: float) -> float:
        fraction = (
            (self.y_zero_px - float(y_px))
            / (self.y_zero_px - self.y_reference_px)
        )
        return self.y_zero_value + fraction * (
            self.y_reference_value - self.y_zero_value
        )


AXES = {
    # Figure 6.5 has three vertically tiled panels with identical scales.
    "6.5a": Axis(861.0, 1400.0, 0.0, 13.0, 829.0, 260.0, 0.0, 7.5),
    "6.5b": Axis(861.0, 1400.0, 0.0, 13.0, 1398.0, 829.0, 0.0, 7.5),
    "6.5c": Axis(861.0, 1400.0, 0.0, 13.0, 1967.0, 1398.0, 0.0, 7.5),
    "6.6": Axis(751.0, 1555.0, 0.0, 65.0, 1108.0, 257.0, 0.0, 9.0),
    "6.9a": Axis(904.0, 1529.0, 0.0, 0.005, 946.0, 284.0, 0.0, 0.005),
    "6.9b": Axis(905.0, 1530.0, 0.0, 0.005, 1827.0, 1166.0, 0.0, 0.005),
    "6.14a": Axis(820.0, 1474.0, -5.0, 45.0, 919.0, 268.0, 0.0, 1.1),
    "6.14b": Axis(820.0, 1474.0, -5.0, 45.0, 1770.0, 1120.0, 0.0, 20.0),
    "6.15a": Axis(831.0, 1467.0, -5.0, 45.0, 928.0, 260.0, 0.0, 1.1),
    "6.15b": Axis(831.0, 1467.0, -5.0, 45.0, 1781.0, 1110.0, 0.0, 11.0),
}


@dataclass(frozen=True)
class DepthPoint:
    panel: str
    energy_eV: int
    substrate: str
    deposited_fc_thickness_A: float
    etch_step_s: float
    marker_x_px: float
    marker_y_px: float
    marker_fill: str
    replicate: int = 1

    @property
    def etch_depth_A_per_cycle(self) -> float:
        return AXES[self.panel].y_value(self.marker_y_px)


def _fig65(panel, energy, substrate, points):
    return tuple(
        DepthPoint(
            panel=panel,
            energy_eV=energy,
            substrate=substrate,
            deposited_fc_thickness_A=AXES[panel].x_value(x_px),
            etch_step_s=40.0,
            marker_x_px=x_px,
            marker_y_px=y_px,
            marker_fill="filled",
            replicate=replicate,
        )
        for x_px, y_px, replicate in points
    )


# Marker centers from the 240-dpi source raster.  Red circles were isolated
# with a full-resolution RGB mask and a distance transform; filled black
# squares were isolated independently from the grayscale mask.  Connected
# overlapping markers in Fig. 6.5b were retained as separate local maxima and
# reconciled against the source page at original resolution.
DEPTH_POINTS = (
    *_fig65("6.5a", 20, "Si", ((1018.2, 662.0, 1),)),
    *_fig65(
        "6.5a",
        20,
        "SiO2",
        (
            (1018.0, 829.0, 1),
            (1101.0, 733.5, 1),
            (1120.5, 775.5, 1),
            (1215.5, 637.0, 1),
        ),
    ),
    *_fig65(
        "6.5b",
        25,
        "Si",
        (
            (989.2, 1215.5, 1),
            (1048.4, 1156.8, 1),
            (1061.4, 1121.4, 1),
            (1061.4, 1217.6, 2),
            (1159.2, 1033.9, 1),
            (1335.4, 1015.7, 1),
        ),
    ),
    *_fig65(
        "6.5b",
        25,
        "SiO2",
        (
            (978.5, 1317.7, 1),
            (996.9, 1273.0, 1),
            (1040.3, 1239.2, 1),
            (1060.0, 1287.0, 1),
            (1077.0, 1283.0, 1),
            (1083.0, 1255.0, 1),
            (1154.5, 1213.0, 1),
            (1196.0, 1233.0, 1),
        ),
    ),
    *_fig65(
        "6.5c",
        30,
        "Si",
        (
            (1010.1, 1724.0, 1),
            (1078.1, 1726.7, 1),
            (1078.3, 1656.0, 2),
            (1093.8, 1699.8, 1),
            (1126.2, 1443.5, 1),
        ),
    ),
    *_fig65(
        "6.5c",
        30,
        "SiO2",
        (
            (977.9, 1872.0, 1),
            (1024.9, 1745.1, 1),
            (1048.0, 1621.5, 1),
            (1139.0, 1645.0, 1),
            (1265.0, 1628.5, 1),
        ),
    ),
    # Figure 6.6: the x values are the printed 20/40/60 s conditions.
    DepthPoint("6.6", 25, "SiO2", 5.0, 20.0, 998.0, 1029.0, "filled"),
    DepthPoint("6.6", 25, "SiO2", 5.0, 40.0, 1245.0, 925.0, "filled"),
    DepthPoint("6.6", 25, "SiO2", 5.0, 60.0, 1492.0, 841.0, "filled"),
    DepthPoint("6.6", 25, "Si", 5.0, 20.0, 997.5, 913.5, "filled"),
    DepthPoint("6.6", 25, "Si", 5.0, 40.0, 1245.5, 863.5, "filled"),
    DepthPoint("6.6", 25, "Si", 5.0, 60.0, 1491.5, 1009.5, "filled"),
    DepthPoint("6.6", 30, "SiO2", 5.0, 20.0, 999.0, 856.0, "open"),
    DepthPoint("6.6", 30, "SiO2", 5.0, 40.0, 1245.0, 679.0, "open"),
    DepthPoint("6.6", 30, "SiO2", 5.0, 60.0, 1493.0, 726.5, "open"),
    DepthPoint("6.6", 30, "Si", 5.0, 20.0, 999.5, 815.0, "open"),
    DepthPoint("6.6", 30, "Si", 5.0, 40.0, 1245.5, 772.5, "open", 1),
    # The source raster contains a second open Si marker at the same 40 s
    # condition.  It is retained as a plotted replicate, not averaged away.
    DepthPoint("6.6", 30, "Si", 5.0, 40.0, 1245.5, 809.0, "open", 2),
    DepthPoint("6.6", 30, "Si", 5.0, 60.0, 1492.5, 913.0, "open"),
)


@dataclass(frozen=True)
class YieldPoint:
    panel: str
    energy_eV: int
    substrate: str
    etch_step_s: float
    marker_x_px: float
    marker_y_px: float
    marker_fill: str
    replicate: int = 1

    @property
    def fluorine_per_incident_ion(self) -> float:
        return AXES[self.panel].x_value(self.marker_x_px)

    @property
    def substrate_units_per_incident_ion(self) -> float:
        return AXES[self.panel].y_value(self.marker_y_px)


# Figure 6.9 is the source's own cycle-normalized response: F in the
# deposited film per incident Ar+ versus SiO2 formula units or Si atoms
# removed per incident Ar+.  The marker centers below were isolated with PIL
# RGB/grayscale masks and visually reconciled against the 240-dpi page.  Panel
# (a) varies deposited FC thickness at fixed 40 s; the individual thickness
# assignment is intentionally not reverse-mapped from Figure 6.5.  Panel (b)
# varies the printed 20/40/60 s step at fixed 5 A deposition; its duplicate
# open Si marker at 40 s is retained.
YIELD_POINTS = (
    # 6.9a: 25 eV, 40 s, 3--11 A deposited FC thickness sweep.
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1015.5, 903.5, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1032.5, 881.0, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1073.5, 863.5, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1093.0, 888.5, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1115.5, 873.5, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1182.5, 849.0, "filled"),
    YieldPoint("6.9a", 25, "SiO2", 40.0, 1222.0, 859.5, "filled"),
    YieldPoint("6.9a", 25, "Si", 40.0, 1025.5, 771.5, "filled"),
    YieldPoint("6.9a", 25, "Si", 40.0, 1082.0, 714.0, "filled"),
    YieldPoint("6.9a", 25, "Si", 40.0, 1094.5, 681.0, "filled"),
    YieldPoint("6.9a", 25, "Si", 40.0, 1187.0, 596.0, "filled"),
    YieldPoint("6.9a", 25, "Si", 40.0, 1354.0, 580.0, "filled"),
    # 6.9b: 5 A deposited FC film, 20--60 s step-length sweep.
    YieldPoint("6.9b", 25, "SiO2", 60.0, 1052.5, 1752.5, "filled"),
    YieldPoint("6.9b", 25, "SiO2", 40.0, 1118.0, 1750.5, "filled"),
    YieldPoint("6.9b", 25, "SiO2", 20.0, 1370.0, 1761.5, "filled"),
    YieldPoint("6.9b", 30, "SiO2", 60.0, 1018.0, 1716.0, "open"),
    YieldPoint("6.9b", 30, "SiO2", 40.0, 1083.5, 1641.5, "open"),
    YieldPoint("6.9b", 30, "SiO2", 20.0, 1297.0, 1610.0, "open"),
    YieldPoint("6.9b", 25, "Si", 60.0, 1035.0, 1781.0, "filled"),
    YieldPoint("6.9b", 25, "Si", 40.0, 1111.0, 1631.0, "filled"),
    YieldPoint("6.9b", 25, "Si", 20.0, 1271.0, 1529.0, "filled"),
    YieldPoint("6.9b", 30, "Si", 60.0, 1016.5, 1729.5, "open"),
    YieldPoint("6.9b", 30, "Si", 40.0, 1111.0, 1593.5, "open", 1),
    YieldPoint("6.9b", 30, "Si", 40.0, 1123.5, 1571.0, "open", 2),
    YieldPoint("6.9b", 30, "Si", 20.0, 1432.0, 1377.5, "open"),
)


@dataclass(frozen=True)
class XPSPoint:
    panel: str
    film_thickness_A: int
    observable: str
    time_s: float
    marker_x_px: float
    marker_y_px: float

    @property
    def value(self) -> float:
        return AXES[self.panel].y_value(self.marker_y_px)


def _xps(panel, thickness, observable, points):
    return tuple(
        XPSPoint(panel, thickness, observable, time, x_px, y_px)
        for time, x_px, y_px in points
    )


XPS_POINTS = (
    *_xps(
        "6.14a", 5, "film_F_over_C_from_C1s",
        ((0, 887, 667), (5, 954, 826), (15, 1083, 819), (40, 1408, 838)),
    ),
    *_xps(
        "6.14a", 5, "delta_F_over_C_substrate_proxy",
        ((0, 887, 697), (5, 954, 718), (15, 1083, 432), (40, 1408, 283)),
    ),
    *_xps(
        "6.14b", 5, "CFx_C1s_intensity_kcps",
        ((0, 886, 1470), (5, 951, 1670), (15, 1081, 1700), (40, 1407, 1708)),
    ),
    *_xps(
        "6.14b", 5, "F1s_intensity_kcps",
        ((0, 886, 1204), (5, 951, 1454), (15, 1080, 1367), (40, 1407, 1227)),
    ),
    *_xps(
        "6.15a", 11, "film_F_over_C_from_C1s",
        ((0, 895, 622), (5, 959, 712), (15, 1084, 756), (40, 1400, 808)),
    ),
    *_xps(
        "6.15a", 11, "delta_F_over_C_substrate_proxy",
        ((0, 895, 909), (5, 959, 891), (15, 1085, 911), (40, 1400, 901)),
    ),
    *_xps(
        "6.15b", 11, "CFx_C1s_intensity_kcps",
        ((0, 893, 1186), (5, 959, 1358), (15, 1085, 1433), (40, 1404, 1546)),
    ),
    *_xps(
        "6.15b", 11, "F1s_intensity_kcps",
        ((0, 895, 1149), (5, 959, 1287), (15, 1085, 1399), (40, 1404, 1495)),
    ),
)


DEPTH_FIELDS = (
    "source_figure", "panel", "energy_eV", "substrate",
    "deposited_fc_thickness_A", "etch_step_s", "etch_depth_A_per_cycle",
    "replicate", "marker_fill", "marker_center_x_px", "marker_center_y_px",
    "digitization_uncertainty_depth_A", "experimental_uncertainty_reported",
    "film_thickness_semantics", "boundary_evidence_tier", "observation_id",
    "source_pdf_sha256",
)
XPS_FIELDS = (
    "source_figure", "panel", "film_thickness_A", "energy_eV",
    "etch_step_s", "time_s", "observable", "value", "units",
    "marker_center_x_px", "marker_center_y_px", "digitization_uncertainty",
    "experimental_uncertainty_reported", "quantity_semantics",
    "supports_absolute_areal_atom_inventory", "boundary_evidence_tier",
    "observation_id", "source_pdf_sha256",
)
YIELD_FIELDS = (
    "source_figure", "panel", "energy_eV", "substrate",
    "sweep_variable", "deposited_fc_thickness_A", "etch_step_s",
    "fluorine_per_incident_ion", "substrate_units_per_incident_ion",
    "substrate_unit_semantics", "replicate", "marker_fill",
    "marker_center_x_px", "marker_center_y_px",
    "digitization_uncertainty_ratio", "experimental_uncertainty_reported",
    "source_derived_normalization", "supports_wafer_flux_reconstruction",
    "boundary_evidence_tier", "observation_id", "source_pdf_sha256",
)


def depth_rows():
    rows = []
    for index, point in enumerate(DEPTH_POINTS, 1):
        if point.panel == "6.6":
            thickness = point.deposited_fc_thickness_A
        else:
            thickness = AXES[point.panel].x_value(point.marker_x_px)
        rows.append({
            "source_figure": f"Figure {point.panel.split('a')[0].split('b')[0].split('c')[0]}",
            "panel": point.panel,
            "energy_eV": str(point.energy_eV),
            "substrate": point.substrate,
            "deposited_fc_thickness_A": f"{thickness:.4f}",
            "etch_step_s": f"{point.etch_step_s:.1f}",
            "etch_depth_A_per_cycle": f"{point.etch_depth_A_per_cycle:.4f}",
            "replicate": str(point.replicate),
            "marker_fill": point.marker_fill,
            "marker_center_x_px": f"{point.marker_x_px:.1f}",
            "marker_center_y_px": f"{point.marker_y_px:.1f}",
            "digitization_uncertainty_depth_A": "0.05",
            "experimental_uncertainty_reported": "false",
            "film_thickness_semantics": (
                "ellipsometric_optical_thickness_not_absolute_C_F_inventory"
            ),
            "boundary_evidence_tier": "A_direct_surface_measurement",
            "observation_id": f"metzler_depth_{index:03d}",
            "source_pdf_sha256": SOURCE["pdf_sha256"],
        })
    return rows


def xps_rows():
    rows = []
    for index, point in enumerate(XPS_POINTS, 1):
        ratio = "over_C" in point.observable
        rows.append({
            "source_figure": f"Figure {point.panel[:-1]}",
            "panel": point.panel,
            "film_thickness_A": str(point.film_thickness_A),
            "energy_eV": "25",
            "etch_step_s": "40.0",
            "time_s": f"{point.time_s:.1f}",
            "observable": point.observable,
            "value": f"{point.value:.4f}",
            "units": "ratio" if ratio else "kcps",
            "marker_center_x_px": f"{point.marker_x_px:.1f}",
            "marker_center_y_px": f"{point.marker_y_px:.1f}",
            "digitization_uncertainty": "0.01" if ratio else "0.10",
            "experimental_uncertainty_reported": "false",
            "quantity_semantics": (
                "relative_XPS_ratio_not_absolute_substrate_F_inventory"
                if point.observable == "delta_F_over_C_substrate_proxy"
                else (
                    "XPS_film_composition_ratio_not_absolute_areal_inventory"
                    if ratio
                    else "background_sensitive_XPS_peak_intensity"
                )
            ),
            "supports_absolute_areal_atom_inventory": "false",
            "boundary_evidence_tier": "A_direct_surface_measurement",
            "observation_id": f"metzler_xps_{index:03d}",
            "source_pdf_sha256": SOURCE["pdf_sha256"],
        })
    return rows


def yield_rows():
    rows = []
    for index, point in enumerate(YIELD_POINTS, 1):
        panel_a = point.panel == "6.9a"
        rows.append({
            "source_figure": "Figure 6.9",
            "panel": point.panel,
            "energy_eV": str(point.energy_eV),
            "substrate": point.substrate,
            "sweep_variable": (
                "deposited_fc_thickness_3_to_11_A"
                if panel_a else "etch_step_length_20_to_60_s"
            ),
            "deposited_fc_thickness_A": "" if panel_a else "5.0",
            "etch_step_s": f"{point.etch_step_s:.1f}",
            "fluorine_per_incident_ion": (
                f"{point.fluorine_per_incident_ion:.7f}"
            ),
            "substrate_units_per_incident_ion": (
                f"{point.substrate_units_per_incident_ion:.7f}"
            ),
            "substrate_unit_semantics": (
                "SiO2_formula_units_per_Ar_ion"
                if point.substrate == "SiO2"
                else "Si_atoms_per_Ar_ion"
            ),
            "replicate": str(point.replicate),
            "marker_fill": point.marker_fill,
            "marker_center_x_px": f"{point.marker_x_px:.1f}",
            "marker_center_y_px": f"{point.marker_y_px:.1f}",
            "digitization_uncertainty_ratio": "0.00002",
            "experimental_uncertainty_reported": "false",
            "source_derived_normalization": (
                "author_cycle_average_from_removed_atoms_and_incident_ions;"
                "F_numerator_from_deposited_FC_film"
            ),
            "supports_wafer_flux_reconstruction": "false",
            "boundary_evidence_tier": "A_direct_surface_response_source_derived",
            "observation_id": f"metzler_yield_{index:03d}",
            "source_pdf_sha256": SOURCE["pdf_sha256"],
        })
    return rows


def _csv_text(fields, rows):
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def _sha_text(text):
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _manifest(depth_text, yield_text, xps_text):
    thin_film = [
        point for point in XPS_POINTS
        if point.panel == "6.14a"
        and point.observable == "delta_F_over_C_substrate_proxy"
    ]
    thick_film = [
        point for point in XPS_POINTS
        if point.panel == "6.15a"
        and point.observable == "delta_F_over_C_substrate_proxy"
    ]
    payload = {
        "manifest_id": "METZLER-2016-CYCLIC-FC-ALE-R2",
        "source": SOURCE,
        "digitization": {
            "method": (
                "official UMD thesis PDF; Poppler 240-dpi render; "
                "PIL/NumPy-compatible RGB and grayscale component isolation; "
                "distance-transform marker centers; original-resolution visual "
                "reconciliation and QA overlays"
            ),
            "source_figures_visually_inspected_at_full_resolution": True,
            "retrospective_board_not_value_blind": True,
            "axis_calibration": {
                key: asdict(value) for key, value in AXES.items()
            },
        },
        "experimental_conditions": {
            "source_power_W": 200,
            "pressure_mTorr": 10,
            "argon_flow_sccm": 50,
            "sample_temperature_C": 10,
            "precursor_depletion_before_bias_s": 12,
            "nominal_ion_energy_eV": [20, 25, 30],
            "xps_cycle_condition": {
                "precursor": "C4F8",
                "ion_energy_eV": 25,
                "etch_step_s": 40,
            },
        },
        "claim_boundary": {
            "valid": [
                "cyclic etch depth versus deposited FC optical thickness",
                "cyclic etch depth versus low-energy Ar-ion step duration",
                "cycle-averaged substrate removal versus F per incident Ar ion",
                "time-resolved FC-film F/C and XPS intensity trajectories",
                "relative substrate-fluorination proxy delta-F/C",
                "a finite-thickness film-transfer/mixing model-form test",
            ],
            "not_valid": [
                "absolute C or F areal inventory inferred from optical thickness",
                "absolute substrate F inventory inferred from delta-F/C",
                "wafer flux reconstructed from the cycle-normalized ratios",
                "species-resolved neutral or positive-ion wafer flux",
                "an ion energy-angle distribution",
                "a value-blind or prospectively held-out claim",
                "permission to tune feature depth",
            ],
        },
        "derived_checks": {
            "depth_markers": len(DEPTH_POINTS),
            "yield_markers": len(YIELD_POINTS),
            "xps_markers": len(XPS_POINTS),
            "thin_film_delta_F_over_C_increase_0_to_40s": round(
                thin_film[-1].value - thin_film[0].value, 4
            ),
            "thick_film_delta_F_over_C_max": round(
                max(point.value for point in thick_film), 4
            ),
            "thin_film_transfers_more_F_to_substrate_proxy": (
                thin_film[-1].value
                > 10.0 * max(point.value for point in thick_film)
            ),
        },
        "outputs": {
            "depth_csv": {
                "path": str(DEPTH_CSV.relative_to(ROOT)),
                "sha256": _sha_text(depth_text),
            },
            "yield_csv": {
                "path": str(YIELD_CSV.relative_to(ROOT)),
                "sha256": _sha_text(yield_text),
            },
            "xps_csv": {
                "path": str(XPS_CSV.relative_to(ROOT)),
                "sha256": _sha_text(xps_text),
            },
        },
    }
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def expected_files():
    depth_text = _csv_text(DEPTH_FIELDS, depth_rows())
    yield_text = _csv_text(YIELD_FIELDS, yield_rows())
    xps_text = _csv_text(XPS_FIELDS, xps_rows())
    return {
        DEPTH_CSV: depth_text,
        YIELD_CSV: yield_text,
        XPS_CSV: xps_text,
        MANIFEST_PATH: _manifest(depth_text, yield_text, xps_text),
    }


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _render_pages(pdf_path, output_dir):
    rendered = {}
    for page in SOURCE["pages"]:
        prefix = output_dir / f"p{page}"
        subprocess.run(
            [
                "pdftoppm", "-f", str(page), "-l", str(page),
                "-r", "240", "-singlefile", "-png",
                str(pdf_path), str(prefix),
            ],
            check=True,
        )
        rendered[page] = prefix.with_suffix(".png")
    return rendered


def _verify_rasters(rendered):
    for page, path in rendered.items():
        image = Image.open(path)
        if image.size != (2040, 2640):
            raise ValueError(
                f"page {page}: expected 2040x2640, got {image.size}"
            )
        expected = SOURCE["pages"][page]["png_sha256"]
        actual = _sha256(path)
        if actual != expected:
            raise ValueError(
                f"page {page}: raster checksum mismatch: {actual}"
            )


def _draw_overlays(rendered, overlay_dir):
    overlay_dir.mkdir(parents=True, exist_ok=True)
    points_by_page = {
        141: [
            point for point in DEPTH_POINTS if point.panel.startswith("6.5")
        ],
        142: [point for point in DEPTH_POINTS if point.panel == "6.6"],
        147: list(YIELD_POINTS),
        157: [point for point in XPS_POINTS if point.panel.startswith("6.14")],
        158: [point for point in XPS_POINTS if point.panel.startswith("6.15")],
    }
    for page, points in points_by_page.items():
        image = Image.open(rendered[page]).convert("RGB")
        draw = ImageDraw.Draw(image)
        for point in points:
            x = point.marker_x_px
            y = point.marker_y_px
            radius = 8
            draw.ellipse(
                (x - radius, y - radius, x + radius, y + radius),
                outline=(0, 120, 255),
                width=3,
            )
            draw.line((x - 12, y, x + 12, y), fill=(0, 120, 255), width=1)
            draw.line((x, y - 12, x, y + 12), fill=(0, 120, 255), width=1)
        image.save(overlay_dir / f"p{page}_qa.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    parser.add_argument("--overlay-dir", type=Path)
    parser.add_argument(
        "--write",
        action="store_true",
        help="write the replayed CSVs/manifest instead of comparing them",
    )
    args = parser.parse_args()

    if not args.pdf.exists():
        raise FileNotFoundError(
            f"download the official thesis PDF to {args.pdf}"
        )
    actual_pdf_sha = _sha256(args.pdf)
    if actual_pdf_sha != SOURCE["pdf_sha256"]:
        raise ValueError(f"source PDF checksum mismatch: {actual_pdf_sha}")

    with tempfile.TemporaryDirectory(prefix="metzler_2016_") as temp:
        rendered = _render_pages(args.pdf, Path(temp))
        _verify_rasters(rendered)
        if args.overlay_dir:
            _draw_overlays(rendered, args.overlay_dir)

    for path, text in expected_files().items():
        if args.write:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text, encoding="utf-8")
        elif path.read_text(encoding="utf-8") != text:
            raise AssertionError(f"committed file differs from replay: {path}")

    print(
        f"verified {len(DEPTH_POINTS)} depth markers, "
        f"{len(YIELD_POINTS)} yield markers, {len(XPS_POINTS)} XPS markers, "
        f"and five full-resolution pages"
    )


if __name__ == "__main__":
    main()

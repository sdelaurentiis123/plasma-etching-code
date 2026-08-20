"""Verified access to the Sayyed et al. Bosch process time series.

The Zenodo record stores each wafer as a NetCDF4/HDF5 group.  Integer samples
index one shared lossless dictionary.  This module keeps that decoding separate
from experimental depth labels: reactor inputs can be loaded and summarized
without opening either wafer-measurement CSV.

``Gas5Flow`` and ``Gas4Flow`` are identified here as the SF6 and C4F8 channels,
respectively.  That mapping is an explicit inference, not metadata printed in
the source file: Gas5 has the 600-sccm, approximately 4.5-second, 100-cycle
etch-phase waveform while Gas4 has the 300-sccm, approximately 1.5-second,
100-cycle passivation waveform documented by the source README.  Every summary
records the inference and the thresholds used to make it auditable.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import md5
import math
from pathlib import Path
import re
from types import MappingProxyType
from typing import Mapping

import numpy as np


PROCESS_DATA_MD5 = "4567d24ec2125102a2e5129203ba31fa"
PROCESS_DICTIONARY_MD5 = "0dde5a3a913eb1fa8512ef2f8748fb34"
WAFER_MEASUREMENT_89_POINT_MD5 = "446e75b040eea37b634eeb8f763a62fc"
WAFER_MEASUREMENT_9_POINT_MD5 = "78515caf25e29e558e1859b92f8a4827"

SF6_FLOW_CHANNEL = "Stat3_Etch_MV_Gas5Flow"
C4F8_FLOW_CHANNEL = "Stat3_Etch_MV_Gas4Flow"
SF6_PHASE_THRESHOLD = 300.0
C4F8_PHASE_THRESHOLD = 150.0

_GROUP = re.compile(r"^Day_(\d{4})_(\d{2})_(\d{2})_Wafer_(\d{2})$")
_REQUIRED_CHANNELS = (
    "Stat3_Etch_MV_EpdIntensity",
    "Stat3_Etch_MV_ForeLinePressure",
    "Stat3_Etch_MV_Gas1Flow",
    "Stat3_Etch_MV_Gas2Flow",
    "Stat3_Etch_MV_Gas3Flow",
    C4F8_FLOW_CHANNEL,
    SF6_FLOW_CHANNEL,
    "Stat3_Etch_MV_Gas7Flow",
    "Stat3_Etch_MV_Gas8Flow",
    "Stat3_Etch_MV_Heater1Temp",
    "Stat3_Etch_MV_Heater2Temp",
    "Stat3_Etch_MV_Heater3Temp",
    "Stat3_Etch_MV_Heater4Temp",
    "Stat3_Etch_MV_HeliumBPFlow",
    "Stat3_Etch_MV_HeliumBPPressure",
    "Stat3_Etch_MV_PlatenDcBias",
    "Stat3_Etch_MV_PlatenRFLoadCapacitor",
    "Stat3_Etch_MV_PlatenRFLoadPower",
    "Stat3_Etch_MV_PlatenRFPeakToPeak",
    "Stat3_Etch_MV_PlatenRFReflectedPower",
    "Stat3_Etch_MV_PlatenRFTuningCapacitor",
    "Stat3_Etch_MV_Pressure",
    "Stat3_Etch_MV_SourceRF2LoadPower",
    "Stat3_Etch_MV_SourceRF2PeakToPeak",
    "Stat3_Etch_MV_SourceRF2ReflectedPower",
    "Stat3_Etch_MV_SourceRF2TuningCapacitor",
    "Stat3_Etch_MV_SourceRFLoadPower",
    "Stat3_Etch_MV_SourceRFPeakToPeak",
    "Stat3_Etch_MV_SourceRFReflectedPower",
    "Stat3_Etch_MV_SourceRFTuningCapacitor",
    "Stat3_Etch_MV_moriInnerCurrent",
)

_PHASE_DIAGNOSTICS = MappingProxyType({
    "pressure": "Stat3_Etch_MV_Pressure",
    "source_load_power": "Stat3_Etch_MV_SourceRFLoadPower",
    "source_reflected_power": "Stat3_Etch_MV_SourceRFReflectedPower",
    "platen_load_power": "Stat3_Etch_MV_PlatenRFLoadPower",
    "platen_reflected_power": "Stat3_Etch_MV_PlatenRFReflectedPower",
    "platen_peak_to_peak": "Stat3_Etch_MV_PlatenRFPeakToPeak",
    "platen_dc_bias": "Stat3_Etch_MV_PlatenDcBias",
    "helium_backside_pressure": "Stat3_Etch_MV_HeliumBPPressure",
    "source_current": "Stat3_Etch_MV_moriInnerCurrent",
    "epd_intensity": "Stat3_Etch_MV_EpdIntensity",
})


def _checksum(path: Path) -> str:
    digest = md5()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _readonly(values, *, dtype=float):
    output = np.asarray(values, dtype=dtype).copy()
    output.setflags(write=False)
    return output


@dataclass(frozen=True)
class BoschProcessTrace:
    """One decoded process record, intentionally containing no etch outcome."""

    experiment_key: str
    source_group: str
    process_date: str
    wafer_number: int
    elapsed_s: np.ndarray
    channels: Mapping[str, np.ndarray]

    def __post_init__(self):
        elapsed = _readonly(self.elapsed_s)
        channels = {
            str(name): _readonly(values)
            for name, values in dict(self.channels).items()
        }
        if (
            not self.experiment_key
            or not self.source_group
            or not self.process_date
            or self.wafer_number <= 0
            or elapsed.ndim != 1
            or elapsed.size < 3000
            or elapsed[0] != 0.0
            or np.any(~np.isfinite(elapsed))
            or np.any(np.diff(elapsed) <= 0.0)
            or set(_REQUIRED_CHANNELS) - set(channels)
            or any(value.shape != elapsed.shape for value in channels.values())
            or any(np.any(~np.isfinite(value)) for value in channels.values())
        ):
            raise ValueError("invalid Bosch process trace")
        sample_period = float(np.median(np.diff(elapsed)))
        if not 0.19 <= sample_period <= 0.21:
            raise ValueError("Bosch process record is not sampled at the declared 5 Hz")
        object.__setattr__(self, "elapsed_s", elapsed)
        object.__setattr__(self, "channels", MappingProxyType(channels))


@dataclass(frozen=True)
class BoschProcessSummary:
    """Deterministic, label-free reactor-input features for one wafer."""

    experiment_key: str
    source_group: str
    process_date: str
    wafer_number: int
    metrics: Mapping[str, float | int]

    def __post_init__(self):
        metrics = dict(self.metrics)
        if not metrics or any(
            not np.isfinite(float(value)) for value in metrics.values()
        ):
            raise ValueError("invalid Bosch process summary")
        object.__setattr__(self, "metrics", MappingProxyType(metrics))


@dataclass(frozen=True)
class BoschWaferMeasurementMap:
    """One identified wafer's spatial Si/oxide measurements in micrometres."""

    experiment_key: str
    lot_number: int
    wafer_number: int
    x_um: np.ndarray
    y_um: np.ndarray
    preoxide_thickness_um: np.ndarray
    postoxide_thickness_um: np.ndarray
    stepheight_um: np.ndarray
    oxide_loss_um: np.ndarray
    silicon_depth_um: np.ndarray

    def __post_init__(self):
        arrays = {
            name: _readonly(getattr(self, name))
            for name in (
                "x_um", "y_um", "preoxide_thickness_um",
                "postoxide_thickness_um", "stepheight_um", "oxide_loss_um",
                "silicon_depth_um",
            )
        }
        size = arrays["x_um"].size
        expected_suffix = f"_{self.wafer_number:02d}"
        if (
            not self.experiment_key.endswith(expected_suffix)
            or self.lot_number <= 0
            or self.wafer_number <= 0
            or size == 0
            or any(value.ndim != 1 or value.size != size for value in arrays.values())
            or any(np.any(~np.isfinite(value)) for value in arrays.values())
            or np.any(arrays["preoxide_thickness_um"] <= 0.0)
            or np.any(arrays["postoxide_thickness_um"] < 0.0)
            or np.any(arrays["stepheight_um"] <= 0.0)
            or np.any(arrays["oxide_loss_um"] < 0.0)
            or np.any(arrays["silicon_depth_um"] <= 0.0)
            or len(set(zip(arrays["x_um"], arrays["y_um"]))) != size
        ):
            raise ValueError("invalid Bosch wafer measurement map")
        for name, value in arrays.items():
            object.__setattr__(self, name, value)

    @property
    def wafer_mean_silicon_depth_um(self):
        return float(np.mean(self.silicon_depth_um))

    @property
    def wafer_mean_oxide_loss_um(self):
        return float(np.mean(self.oxide_loss_um))

    @property
    def wafer_mean_selectivity(self):
        oxide = self.wafer_mean_oxide_loss_um
        return self.wafer_mean_silicon_depth_um / oxide if oxide > 0.0 else math.inf


_MEASUREMENT_89_COLUMNS = (
    "experiment_key", "lot_number", "wafer_number", "X", "Y",
    "preox_thickness", "postox_thickness", "postox_thickness_nan",
    "stepheight", "oxide_etch", "si_etch",
)


def load_bosch_wafer_measurements_89pt(
        path, *, allowed_experiment_keys) -> tuple[BoschWaferMeasurementMap, ...]:
    """Load a calibration-only or revealed 89-point measurement asset.

    ``allowed_experiment_keys`` is mandatory.  A row outside that allowlist is
    rejected immediately after reading its string key and before converting
    any numeric outcome field.  Fit code must therefore consume the extracted
    calibration-only asset, never the mixed official source CSV.
    """
    path = Path(path)
    allowed = frozenset(str(key) for key in allowed_experiment_keys)
    if not allowed or any(not key for key in allowed):
        raise ValueError("measurement loader requires a nonempty key allowlist")

    grouped: dict[str, list[dict[str, str]]] = {}
    identities: dict[str, tuple[int, int]] = {}
    with path.open("r", newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != _MEASUREMENT_89_COLUMNS:
            raise ValueError("unexpected Bosch 89-point measurement schema")
        for row_number, row in enumerate(reader, start=2):
            key = str(row["experiment_key"]).strip()
            if key not in allowed:
                raise ValueError(
                    f"measurement row {row_number} is outside the allowed key set: {key!r}")
            try:
                lot = int(row["lot_number"])
                wafer = int(row["wafer_number"])
                numeric = {
                    name: float(row[column])
                    for name, column in (
                        ("x_um", "X"), ("y_um", "Y"),
                        ("preoxide_thickness_um", "preox_thickness"),
                        ("postoxide_thickness_um", "postox_thickness"),
                        ("stepheight_um", "stepheight"),
                        ("oxide_loss_um", "oxide_etch"),
                        ("silicon_depth_um", "si_etch"),
                    )
                }
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"invalid numeric Bosch measurement at row {row_number}") from error
            identity = (lot, wafer)
            if key in identities and identities[key] != identity:
                raise ValueError(f"inconsistent Bosch wafer identity: {key}")
            identities[key] = identity
            grouped.setdefault(key, []).append(numeric)

    output = []
    for key in sorted(grouped):
        rows = grouped[key]
        lot, wafer = identities[key]
        if len(rows) != 89:
            raise ValueError(f"Bosch 89-point wafer is incomplete: {key}")
        output.append(BoschWaferMeasurementMap(
            experiment_key=key, lot_number=lot, wafer_number=wafer,
            **{
                name: np.asarray([row[name] for row in rows], dtype=float)
                for name in rows[0]
            }))
    if not output:
        raise ValueError("no allowed Bosch wafer measurements were found")
    return tuple(output)


def load_bosch_process_traces(
    process_path, dictionary_path, *, verify_checksum: bool = True,
) -> tuple[BoschProcessTrace, ...]:
    """Decode every source group while never touching wafer outcome files.

    ``h5py`` is an optional data-ingestion dependency because the runtime etch
    engine does not otherwise require HDF5.
    """

    try:
        import h5py
    except ImportError as error:  # pragma: no cover - environment dependent
        raise RuntimeError(
            "Bosch process ingestion requires the optional 'reactor-data' dependency"
        ) from error

    process_path = Path(process_path)
    dictionary_path = Path(dictionary_path)
    if verify_checksum:
        if _checksum(process_path) != PROCESS_DATA_MD5:
            raise ValueError(f"checksum mismatch for Bosch process data: {process_path}")
        if _checksum(dictionary_path) != PROCESS_DICTIONARY_MD5:
            raise ValueError(
                f"checksum mismatch for Bosch process dictionary: {dictionary_path}")

    with h5py.File(dictionary_path, "r") as source:
        if set(source) != {"data", "x"}:
            raise ValueError("unexpected Bosch process dictionary schema")
        decoder = np.asarray(source["data"][:], dtype=float)
    if decoder.ndim != 1 or decoder.size != 49290 or np.any(~np.isfinite(decoder)):
        raise ValueError("invalid Bosch process dictionary")

    traces = []
    with h5py.File(process_path, "r") as source:
        for group_name in sorted(source):
            match = _GROUP.fullmatch(group_name)
            if match is None:
                raise ValueError(f"unexpected Bosch process group: {group_name}")
            year, month, day, wafer = match.groups()
            process_date = f"{year}-{month}-{day}"
            wafer_number = int(wafer)
            group = source[group_name]
            if not {"data", "feature", "times"}.issubset(group):
                raise ValueError(f"incomplete Bosch process group: {group_name}")
            names = tuple(
                value.decode("utf-8") if isinstance(value, bytes) else str(value)
                for value in group["feature"][:]
            )
            if len(names) != len(set(names)) or set(_REQUIRED_CHANNELS) - set(names):
                raise ValueError(f"unexpected Bosch feature schema: {group_name}")
            encoded = np.asarray(group["data"][:])
            timestamps = np.asarray(group["times"][:], dtype=float)
            if (
                encoded.ndim != 2
                or encoded.shape != (timestamps.size, len(names))
                or encoded.size == 0
                or int(encoded.max()) >= decoder.size
            ):
                raise ValueError(f"invalid Bosch encoded samples: {group_name}")
            decoded = decoder[encoded]
            traces.append(BoschProcessTrace(
                experiment_key=f"{process_date}_{wafer_number:02d}",
                source_group=group_name,
                process_date=process_date,
                wafer_number=wafer_number,
                elapsed_s=timestamps - timestamps[0],
                channels={name: decoded[:, index] for index, name in enumerate(names)},
            ))
    if len(traces) != 96 or len({trace.experiment_key for trace in traces}) != 96:
        raise ValueError("Bosch process file does not contain the expected 96 wafer records")
    return tuple(traces)


def _episode_count(active: np.ndarray) -> int:
    padded = np.concatenate(([False], np.asarray(active, dtype=bool), [False]))
    return int(np.count_nonzero(np.diff(padded.astype(np.int8)) == 1))


def _phase_statistics(values: np.ndarray, active: np.ndarray):
    selected = np.asarray(values, dtype=float)[active]
    if selected.size == 0:
        raise ValueError("empty Bosch process phase")
    return {
        "mean": float(np.mean(selected)),
        "std": float(np.std(selected)),
        "q10": float(np.quantile(selected, 0.10)),
        "q50": float(np.quantile(selected, 0.50)),
        "q90": float(np.quantile(selected, 0.90)),
        "rms": float(np.sqrt(np.mean(selected * selected))),
    }


def summarize_bosch_process_trace(trace: BoschProcessTrace) -> BoschProcessSummary:
    """Reduce one trace to a stable, units-preserving reactor feature row."""

    if not isinstance(trace, BoschProcessTrace):
        raise TypeError("trace must be a BoschProcessTrace")
    sample_period = float(np.median(np.diff(trace.elapsed_s)))
    sf6 = trace.channels[SF6_FLOW_CHANNEL]
    c4f8 = trace.channels[C4F8_FLOW_CHANNEL]
    sf6_phase = sf6 > SF6_PHASE_THRESHOLD
    c4f8_phase = c4f8 > C4F8_PHASE_THRESHOLD
    if np.any(sf6_phase & c4f8_phase):
        raise ValueError("Bosch etch and passivation phases overlap above thresholds")

    metrics: dict[str, float | int] = {
        "sample_count": int(trace.elapsed_s.size),
        "record_duration_s": float(trace.elapsed_s[-1]),
        "median_sample_period_s": sample_period,
        "sf6_episode_count": _episode_count(sf6_phase),
        "c4f8_episode_count": _episode_count(c4f8_phase),
        "sf6_above_threshold_s": float(np.count_nonzero(sf6_phase) * sample_period),
        "c4f8_above_threshold_s": float(np.count_nonzero(c4f8_phase) * sample_period),
        "transition_or_idle_s": float(
            np.count_nonzero(~(sf6_phase | c4f8_phase)) * sample_period),
        "sf6_flow_dose_machine_units_s": float(np.trapz(sf6, trace.elapsed_s)),
        "c4f8_flow_dose_machine_units_s": float(np.trapz(c4f8, trace.elapsed_s)),
    }
    if not 100 <= metrics["sf6_episode_count"] <= 102:
        raise ValueError("unexpected number of SF6/ignition episodes")
    if metrics["c4f8_episode_count"] != 100:
        raise ValueError("unexpected number of C4F8 passivation episodes")

    for phase_name, active in (("sf6", sf6_phase), ("c4f8", c4f8_phase)):
        for diagnostic_name, channel_name in _PHASE_DIAGNOSTICS.items():
            for statistic, value in _phase_statistics(
                trace.channels[channel_name], active
            ).items():
                metrics[f"{phase_name}_{diagnostic_name}_{statistic}"] = value

    return BoschProcessSummary(
        experiment_key=trace.experiment_key,
        source_group=trace.source_group,
        process_date=trace.process_date,
        wafer_number=trace.wafer_number,
        metrics=metrics,
    )


def summarize_bosch_process_traces(
    traces: tuple[BoschProcessTrace, ...],
) -> tuple[BoschProcessSummary, ...]:
    summaries = tuple(summarize_bosch_process_trace(trace) for trace in traces)
    if len(summaries) != len({summary.experiment_key for summary in summaries}):
        raise ValueError("duplicate Bosch process summary")
    return summaries


def process_ingestion_manifest():
    return {
        "schema": "petch-zenodo-17122442-process-ingestion-v1",
        "source_record": "https://zenodo.org/records/17122442",
        "process_data_md5": PROCESS_DATA_MD5,
        "process_dictionary_md5": PROCESS_DICTIONARY_MD5,
        "source_sampling_hz": 5.0,
        "expected_wafer_records": 96,
        "sf6_channel": SF6_FLOW_CHANNEL,
        "sf6_phase_threshold": SF6_PHASE_THRESHOLD,
        "c4f8_channel": C4F8_FLOW_CHANNEL,
        "c4f8_phase_threshold": C4F8_PHASE_THRESHOLD,
        "gas_channel_mapping_evidence": (
            "waveform-role inference: Gas5 is the 600-unit, approximately 4.5 s etch "
            "waveform; Gas4 is the 300-unit, approximately 1.5 s passivation waveform"
        ),
        "experimental_outcomes_read": False,
    }

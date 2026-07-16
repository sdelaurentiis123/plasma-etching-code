"""First-principles reduced RF sheath models for plasma-to-feature boundary states."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

EPS0 = 8.8541878128e-12
ECHARGE = 1.602176634e-19
AMU = 1.66053906660e-27


def _collisionless_ion_impact_energies(
        *, voltage, period_s, maximum_frequency_hz, mean_voltage_v,
        ion_mass_amu, Te_eV, thickness_m, phases,
        steps_per_period=512, steps_per_transit=512, max_periods=100.0):
    """Integrate collisionless Child-profile ion trajectories for a periodic voltage law."""
    phases = np.asarray(phases, dtype=float)
    shape = phases.shape
    phase = phases.ravel()
    mass = float(ion_mass_amu) * AMU
    s = float(thickness_m)
    period = float(period_s)
    fastest_period = 1.0 / float(maximum_frequency_hz)
    v0 = bohm_speed(Te_eV, ion_mass_amu)
    # Sample the supplied periodic law rather than assuming a sinusoidal peak.
    probe_time = period * np.arange(2048, dtype=float) / 2048.0
    peak_voltage = max(float(np.max(voltage(probe_time))), float(mean_voltage_v), 0.0)
    vmax = np.sqrt(v0 * v0 + 2.0 * ECHARGE * peak_voltage / mass)
    transit_est = 2.0 * s / max(v0 + vmax, 1e-30)
    dt = min(
        fastest_period / int(steps_per_period),
        transit_est / int(steps_per_transit),
    )
    max_steps = int(np.ceil(float(max_periods) * period / dt))
    entry_time = np.mod(phase, 2.0 * np.pi) * period / (2.0 * np.pi)
    x = np.zeros(phase.size)
    v = np.full(phase.size, v0)
    t = entry_time.copy()
    active = np.ones(phase.size, dtype=bool)
    energy = np.full(phase.size, np.nan)

    def acceleration(xx, tt):
        sheath_voltage = np.maximum(voltage(tt), 0.0)
        xi = np.clip(xx / s, 0.0, 1.0)
        field = (4.0 / 3.0) * sheath_voltage / s * np.cbrt(xi)
        return ECHARGE * field / mass

    for _ in range(max_steps):
        if not active.any():
            break
        idx = np.where(active)[0]
        xa = x[idx]
        va = v[idx]
        ta = t[idx]
        a0 = acceleration(xa, ta)
        vh = va + 0.5 * a0 * dt
        xn = xa + vh * dt
        tn = ta + dt
        an = acceleration(xn, tn)
        vn = vh + 0.5 * an * dt
        crossed = xn >= s
        if crossed.any():
            # Interpolate velocity to the physical electrode-crossing time.
            frac = np.clip(
                (s - xa[crossed])
                / np.maximum(xn[crossed] - xa[crossed], 1e-30),
                0.0, 1.0)
            vcross = va[crossed] + frac * (vn[crossed] - va[crossed])
            hit = idx[crossed]
            energy[hit] = 0.5 * mass * vcross * vcross / ECHARGE
            active[hit] = False
        keep = ~crossed
        x[idx[keep]] = xn[keep]
        v[idx[keep]] = vn[keep]
        t[idx[keep]] = tn[keep]
    if active.any():
        raise RuntimeError("ion sheath transit did not finish within max_periods")
    return energy.reshape(shape)


def bohm_speed(Te_eV, ion_mass_amu):
    return float(np.sqrt(ECHARGE * float(Te_eV) / (float(ion_mass_amu) * AMU)))


def child_langmuir_sheath_thickness(V_dc, Te_eV, ion_mass_amu, density_m3):
    """Collisionless planar Child-law thickness from Bohm ion current and mean sheath voltage."""
    mass = float(ion_mass_amu) * AMU
    current = ECHARGE * float(density_m3) * bohm_speed(Te_eV, ion_mass_amu)
    coefficient = (4.0 / 9.0) * EPS0 * np.sqrt(2.0 * ECHARGE / mass)
    return float(np.sqrt(coefficient * float(V_dc) ** 1.5 / current))


@dataclass(frozen=True)
class CollisionlessRFSheath:
    """Planar collisionless RF sheath with finite ion transit time.

    The Child potential shape is ``Phi(x,t)=Vs(t)*(x/s)^(4/3)`` from sheath edge ``x=0`` to wafer
    ``x=s``. Ions enter at the Bohm speed. This is a reduced physical model—not an instantaneous IEDF
    prescription—and is parameterized only by measurable plasma/process quantities.
    """
    V_dc: float
    V_rf: float
    frequency_hz: float
    Te_eV: float
    ion_mass_amu: float
    density_m3: float | None = None
    thickness_m: float | None = None

    @property
    def thickness(self):
        if self.thickness_m is not None:
            return float(self.thickness_m)
        if self.density_m3 is None:
            raise ValueError("density_m3 or thickness_m is required")
        return child_langmuir_sheath_thickness(
            self.V_dc, self.Te_eV, self.ion_mass_amu, self.density_m3)

    def voltage(self, time_s):
        return self.V_dc + self.V_rf * np.sin(2.0 * np.pi * self.frequency_hz * time_s)

    def ion_impact_energies(self, phases, steps_per_period=512, steps_per_transit=512,
                            max_periods=100.0):
        """Integrate independent ions from sheath edge to wafer; return impact energies in eV.

        Velocity-Verlet is used with a step resolving both RF period and estimated mean transit time.
        Trajectories are vectorized over entry phase and deterministic.
        """
        return _collisionless_ion_impact_energies(
            voltage=self.voltage,
            period_s=1.0 / self.frequency_hz,
            maximum_frequency_hz=self.frequency_hz,
            mean_voltage_v=self.V_dc,
            ion_mass_amu=self.ion_mass_amu,
            Te_eV=self.Te_eV,
            thickness_m=self.thickness,
            phases=phases,
            steps_per_period=steps_per_period,
            steps_per_transit=steps_per_transit,
            max_periods=max_periods,
        )


@dataclass(frozen=True)
class PeriodicSheathVoltage:
    """Fourier representation of a measured or modeled periodic sheath drop.

    ``harmonic_number`` is relative to ``fundamental_frequency_hz``.  Both sine and cosine
    amplitudes are retained so a measured waveform can be reconstructed without silently replacing
    it by a single sinusoid.  ``evidence_kind`` is claim-bearing provenance, not a numerical option.
    """
    fundamental_frequency_hz: float
    dc_v: float
    harmonic_number: np.ndarray
    sine_v: np.ndarray
    cosine_v: np.ndarray
    source: str
    evidence_kind: str = "assumed"

    def __post_init__(self):
        harmonic = np.asarray(self.harmonic_number, dtype=int)
        sine = np.asarray(self.sine_v, dtype=float)
        cosine = np.asarray(self.cosine_v, dtype=float)
        if (not np.isfinite(self.fundamental_frequency_hz)
                or self.fundamental_frequency_hz <= 0.0
                or not np.isfinite(self.dc_v) or self.dc_v < 0.0
                or harmonic.ndim != 1 or sine.shape != harmonic.shape
                or cosine.shape != harmonic.shape or harmonic.size == 0
                or np.any(harmonic <= 0) or len(np.unique(harmonic)) != harmonic.size
                or np.any(~np.isfinite(sine)) or np.any(~np.isfinite(cosine))
                or not str(self.source).strip()
                or self.evidence_kind not in {
                    "measured_sheath_voltage", "validated_reactor_model", "assumed"}):
            raise ValueError("invalid periodic sheath-voltage waveform")
        order = np.argsort(harmonic)
        harmonic = harmonic[order]
        sine = sine[order]
        cosine = cosine[order]
        harmonic.setflags(write=False)
        sine.setflags(write=False)
        cosine.setflags(write=False)
        object.__setattr__(self, "harmonic_number", harmonic)
        object.__setattr__(self, "sine_v", sine)
        object.__setattr__(self, "cosine_v", cosine)

    @property
    def period_s(self):
        return 1.0 / float(self.fundamental_frequency_hz)

    @property
    def maximum_frequency_hz(self):
        return float(self.fundamental_frequency_hz * self.harmonic_number[-1])

    @property
    def supports_predictive_boundary(self):
        return self.evidence_kind in {
            "measured_sheath_voltage", "validated_reactor_model"}

    def voltage(self, time_s):
        time = np.asarray(time_s, dtype=float)
        angle = (2.0 * np.pi * self.fundamental_frequency_hz
                 * time[..., None] * self.harmonic_number)
        return (self.dc_v
                + np.sum(self.sine_v * np.sin(angle) + self.cosine_v * np.cos(angle), axis=-1))

    @classmethod
    def sinusoidal(cls, *, dc_v, amplitude_v, frequency_hz, phase_rad=0.0,
                   source, evidence_kind="assumed"):
        """Construct a one-harmonic waveform without losing the declared phase."""
        return cls(
            fundamental_frequency_hz=float(frequency_hz), dc_v=float(dc_v),
            harmonic_number=np.array([1]),
            sine_v=np.array([float(amplitude_v) * np.cos(float(phase_rad))]),
            cosine_v=np.array([float(amplitude_v) * np.sin(float(phase_rad))]),
            source=source, evidence_kind=evidence_kind)


@dataclass(frozen=True)
class CollisionlessWaveformSheath:
    """Finite-transit Child-profile sheath driven by an arbitrary periodic waveform."""
    waveform: PeriodicSheathVoltage
    Te_eV: float
    ion_mass_amu: float
    density_m3: float | None = None
    thickness_m: float | None = None

    def __post_init__(self):
        values = np.asarray([self.Te_eV, self.ion_mass_amu], dtype=float)
        if (np.any(~np.isfinite(values)) or np.any(values <= 0.0)
                or (self.density_m3 is None) == (self.thickness_m is None)):
            raise ValueError("require positive sheath inputs and exactly one thickness closure")
        if self.density_m3 is not None and (
                not np.isfinite(self.density_m3) or self.density_m3 <= 0.0):
            raise ValueError("density_m3 must be positive")
        if self.thickness_m is not None and (
                not np.isfinite(self.thickness_m) or self.thickness_m <= 0.0):
            raise ValueError("thickness_m must be positive")

    @property
    def thickness(self):
        if self.thickness_m is not None:
            return float(self.thickness_m)
        return child_langmuir_sheath_thickness(
            self.waveform.dc_v, self.Te_eV, self.ion_mass_amu, self.density_m3)

    def voltage(self, time_s):
        return self.waveform.voltage(time_s)

    def ion_impact_energies(self, phases, steps_per_period=512, steps_per_transit=512,
                            max_periods=100.0):
        return _collisionless_ion_impact_energies(
            voltage=self.voltage,
            period_s=self.waveform.period_s,
            maximum_frequency_hz=self.waveform.maximum_frequency_hz,
            mean_voltage_v=self.waveform.dc_v,
            ion_mass_amu=self.ion_mass_amu,
            Te_eV=self.Te_eV,
            thickness_m=self.thickness,
            phases=phases,
            steps_per_period=steps_per_period,
            steps_per_transit=steps_per_transit,
            max_periods=max_periods,
        )

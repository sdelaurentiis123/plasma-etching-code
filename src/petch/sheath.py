"""First-principles reduced RF sheath models for plasma-to-feature boundary states."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

EPS0 = 8.8541878128e-12
ECHARGE = 1.602176634e-19
AMU = 1.66053906660e-27


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
        phases = np.asarray(phases, dtype=float)
        shape = phases.shape; phase = phases.ravel()
        mass = self.ion_mass_amu * AMU; s = self.thickness
        omega = 2.0 * np.pi * self.frequency_hz
        period = 1.0 / self.frequency_hz
        v0 = bohm_speed(self.Te_eV, self.ion_mass_amu)
        vmax = np.sqrt(v0 * v0 + 2.0 * ECHARGE * (self.V_dc + abs(self.V_rf)) / mass)
        transit_est = 2.0 * s / max(v0 + vmax, 1e-30)
        dt = min(period / int(steps_per_period), transit_est / int(steps_per_transit))
        max_steps = int(np.ceil(max_periods * period / dt))
        entry_time = phase / omega
        x = np.zeros(phase.size); v = np.full(phase.size, v0)
        t = entry_time.copy(); active = np.ones(phase.size, dtype=bool)
        energy = np.full(phase.size, np.nan)

        def acceleration(xx, tt):
            voltage = np.maximum(self.voltage(tt), 0.0)
            xi = np.clip(xx / s, 0.0, 1.0)
            field = (4.0 / 3.0) * voltage / s * np.cbrt(xi)
            return ECHARGE * field / mass

        for _ in range(max_steps):
            if not active.any():
                break
            idx = np.where(active)[0]
            xa = x[idx]; va = v[idx]; ta = t[idx]
            a0 = acceleration(xa, ta)
            vh = va + 0.5 * a0 * dt
            xn = xa + vh * dt
            tn = ta + dt
            an = acceleration(xn, tn)
            vn = vh + 0.5 * an * dt
            crossed = xn >= s
            if crossed.any():
                # Interpolate velocity to the physical electrode-crossing time.
                frac = np.clip((s - xa[crossed]) / np.maximum(xn[crossed] - xa[crossed], 1e-30), 0.0, 1.0)
                vcross = va[crossed] + frac * (vn[crossed] - va[crossed])
                hit = idx[crossed]
                energy[hit] = 0.5 * mass * vcross * vcross / ECHARGE
                active[hit] = False
            keep = ~crossed
            x[idx[keep]] = xn[keep]; v[idx[keep]] = vn[keep]; t[idx[keep]] = tn[keep]
        if active.any():
            raise RuntimeError("ion sheath transit did not finish within max_periods")
        return energy.reshape(shape)

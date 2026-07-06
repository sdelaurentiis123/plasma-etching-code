"""Differentiable ALE forward model (torch) — the moat: reverse-mode gradients of etch-per-cycle
through the whole cyclic site-balance ROM, w.r.t. ion energy AND the calibration parameters.

Same Vella-Graves ROM as ale.py, re-expressed in torch so autograd flows end to end. Enables:
  - dEPC/dE  (process sensitivity, sub-ms)
  - dEPC/dtheta for any physical parameter (K_Cl mixing, gamma_Cl, yield slopes) -> gradient-based
    CALIBRATION against measured EPC, instead of one-parameter-at-a-time DoE.
  - inverse design: solve for the ion energy (or recipe knob) that hits a target EPC by descent.

No open feature-scale etch tool exposes chemistry gradients like this (ViennaPS: no ALE; Kushner:
closed/CPU/non-diff). Coarser dt / fewer cycles than ale.py (autograd graph memory) but same physics;
ale.py stays the high-accuracy forward reference.
"""
import torch

# fixed constants (mirror ale.py)
SIGMA1 = 1.0e15
J_CL2 = 9.8e17
J_AR = 3.7e16
N_SI = 5.0e22
E_TH_MIX = 2.81


def _yields(E, params):
    """Energy-dependent yields + mixed-layer ratio, as torch tensors. params overrides the calibration
    coefficients (defaults = paper values) so any of them can carry a gradient."""
    p = params
    z = torch.zeros_like(E)
    return dict(
        Y_Cl=torch.maximum(z, p["yc_a"] * E + p["yc_b"]),
        Y_SiCl=torch.maximum(z, p["ysc_a"] * E + p["ysc_b"]),
        Y_SiCl2=torch.maximum(z, p["ysc2_a"] * E + p["ysc2_b"]),
        Y_Si=torch.maximum(z, p["ysi_a"] * E + p["ysi_b"]),
        s2_ratio=torch.clamp(0.77 * (torch.sqrt(torch.clamp(E, min=0.0)) - E_TH_MIX ** 0.5), min=1e-6),
    )


def default_params(device="cpu", dtype=torch.float64, requires_grad_keys=()):
    """The paper's calibration coefficients as (optionally grad-tracked) tensors."""
    base = dict(yc_a=3.2e-3, yc_b=-4.25e-2, ysc_a=1e-4, ysc_b=-1.4e-3,
                ysc2_a=2e-3, ysc2_b=-2.9e-2, ysi_a=4.5e-5, ysi_b=-9e-4,
                gamma_cl=0.25, k_cl=0.45)
    out = {}
    for k, v in base.items():
        out[k] = torch.tensor(v, device=device, dtype=dtype,
                              requires_grad=(k in requires_grad_keys))
    return out


def epc_torch(E, params=None, n_cycles=2, t_dose=0.112, t_barr=113.5, dt=0.05, t_active=45.0):
    """Differentiable etch-per-cycle (Angstrom/cycle) at ion energy E (torch scalar tensor).
    Speed: the coverages decay within the first tens of seconds, after which only the constant
    bare-Si sputter leak J_Ar*Y_Si*(1-theta2) continues -- so integrate the transient finely for
    t_active seconds, then add the leak for the (t_barr - t_active) tail analytically. Same physics
    as ale.py, ~40x fewer autograd steps. Cyclic steady state after ~2 cycles."""
    if params is None:
        params = default_params(dtype=E.dtype, device=E.device)
    y = _yields(E, params)
    s1 = SIGMA1
    s2 = y["s2_ratio"] * SIGMA1
    alpha = 2.0 * J_CL2 * params["gamma_cl"] / SIGMA1
    K = params["k_cl"]
    theta1 = torch.zeros((), dtype=E.dtype, device=E.device)
    theta2 = torch.zeros((), dtype=E.dtype, device=E.device)
    n = max(1, int(round(min(t_active, t_barr) / dt)))
    step = min(t_active, t_barr) / n
    last_epc = torch.zeros((), dtype=E.dtype, device=E.device)
    for _ in range(n_cycles):
        theta1 = 1.0 - (1.0 - theta1) * torch.exp(-alpha * t_dose)   # modification (analytic)
        N_si = torch.zeros((), dtype=E.dtype, device=E.device)
        for _ in range(n):                                           # transient (forward Euler)
            r1 = (J_AR / s1) * (y["Y_Cl"] * theta1 + y["Y_SiCl"] * theta1
                                + 2.0 * y["Y_SiCl2"] * theta1 ** 2 + K * theta1)
            r2 = (J_AR / s2) * (y["Y_Cl"] * theta2 + y["Y_SiCl"] * theta2
                                + 2.0 * y["Y_SiCl2"] * theta2 ** 2 - K * theta1)
            tsum = theta1 + theta2
            j_si = J_AR * (y["Y_SiCl"] * tsum + y["Y_SiCl2"] * tsum + y["Y_Si"] * (1.0 - theta2))
            theta1 = torch.clamp(theta1 - step * r1, min=0.0)
            theta2 = torch.clamp(theta2 - step * r2, min=0.0, max=1.0)
            N_si = N_si + step * j_si
        if t_barr > t_active:                                        # constant sputter-leak tail
            N_si = N_si + J_AR * y["Y_Si"] * (1.0 - theta2) * (t_barr - t_active)
        last_epc = N_si / N_SI * 1e8
    return last_epc


def dEPC_dE(E_value, **kw):
    """Reverse-mode dEPC/dE at a given energy. Returns (epc, gradient)."""
    E = torch.tensor(float(E_value), dtype=torch.float64, requires_grad=True)
    epc = epc_torch(E, **kw)
    epc.backward()
    return float(epc.detach()), float(E.grad)


def invert_energy_for_epc(target_epc, E0=18.0, steps=15, lr=None, **kw):
    """Inverse design: solve for the Ar+ energy that yields target_epc via gradient descent (Newton
    step on the scalar). Returns (E_solved, epc_at_solution, history)."""
    E = torch.tensor(float(E0), dtype=torch.float64, requires_grad=True)
    hist = []
    for _ in range(steps):
        if E.grad is not None:
            E.grad.zero_()
        epc = epc_torch(E, **kw)
        loss = (epc - target_epc) ** 2
        loss.backward()
        g = E.grad
        with torch.no_grad():
            # Newton-ish on EPC: dE = -(epc-target)/(dEPC/dE); dEPC/dE = grad(loss)/(2(epc-target))
            denom = g / (2.0 * (epc - target_epc) + 1e-30)
            stepE = (epc - target_epc) / (denom + 1e-30)
            E -= torch.clamp(stepE, -3.0, 3.0)
            E.clamp_(13.0, 40.0)
        hist.append((float(E.detach()), float(epc.detach())))
    return float(E.detach()), float(epc_torch(E.detach())), hist

#!/usr/bin/env python3
"""Differentiable inverse design: optimize a RECIPE parameter so the etch matches a TARGET.

The payoff of building in Warp. We take the real per-face flux from the 3D simulator, express the
etch chemistry as a differentiable Warp kernel, and ask: "what ion energy makes the trench floor
etch at a target rate?" wp.Tape gives d(loss)/d(energy); gradient descent finds the recipe.

This is auto-calibration / inverse design in one runnable demo (CPU): target outcome -> recipe.
"""
import numpy as np
import warp as wp
import petch
from petch import threed as t3

wp.init()
DEV = "cpu"
P = petch.PAR


@wp.kernel
def floor_rate_sum(m_i: wp.array(dtype=float), m_F: wp.array(dtype=float),
                   m_O: wp.array(dtype=float), cos_i: wp.array(dtype=float),
                   is_floor: wp.array(dtype=float), E: wp.array(dtype=float),
                   sumV: wp.array(dtype=float)):
    f = wp.tid()
    e = E[0]
    Yie = 7.0 * wp.max(wp.sqrt(e) - wp.sqrt(15.0), 0.0)
    Ysp = 0.0337 * wp.max(wp.sqrt(e) - wp.sqrt(20.0), 0.0)
    Yp = 3.0 * wp.max(wp.sqrt(e) - wp.sqrt(10.0), 0.0)
    fang = wp.clamp(cos_i[f], 0.0, 1.0)
    Fi = 12.0 * m_i[f] * fang
    Fev = 1800.0 * m_F[f]
    Fp = 100.0 * m_O[f]
    rF = 0.2 * Fev / (Yie * Fi + 1.0e-9)
    rO = 0.3 * Fp / (Yp * Fi + 1.0e-9)
    bare = 1.0 / (1.0 + rF + rO)
    thF = rF * bare
    V = (1.0 / 5.02) * (Yie * Fi * thF + Ysp * Fi * bare)
    if is_floor[f] > 0.5:
        wp.atomic_add(sumV, 0, V)


@wp.kernel
def mean_loss(sumV: wp.array(dtype=float), n: float, target: float, loss: wp.array(dtype=float)):
    loss[0] = (sumV[0] / n - target) * (sumV[0] / n - target)


def floor_rate_np(m_i, m_F, m_O, cos_i, floor, E):
    """Same chemistry in numpy — to set a reachable target and to sanity-check."""
    Yie = 7.0 * max(np.sqrt(E) - np.sqrt(15.0), 0)
    Ysp = 0.0337 * max(np.sqrt(E) - np.sqrt(20.0), 0)
    Yp = 3.0 * max(np.sqrt(E) - np.sqrt(10.0), 0)
    fang = np.clip(cos_i, 0, 1)
    Fi = 12.0 * m_i * fang; Fev = 1800.0 * m_F; Fp = 100.0 * m_O
    rF = 0.2 * Fev / (Yie * Fi + 1e-9); rO = 0.3 * Fp / (Yp * Fi + 1e-9)
    bare = 1.0 / (1 + rF + rO); thF = rF * bare
    V = (1 / 5.02) * (Yie * Fi * thF + Ysp * Fi * bare)
    return V[floor].mean()


def main():
    # 1. real flux from the simulator on an etched trench
    par = dict(P); par['rate_scale'] = 0.3
    geo = t3.run_etch_3d(Lx=10, Ly=4, Lz=14, dx=0.4, trench_width=4, mask_th=2, sub_top=10,
                         t_end=2.0, n_steps=10, par=par, flags=petch.Flags(),
                         n_ion=8000, n_neu=8000, verbose=False)
    verts, faces, cent, areas = t3.extract_mesh_3d(geo['phi'], 0.4)
    mesh = wp.Mesh(points=wp.array(verts, dtype=wp.vec3, device=DEV),
                   indices=wp.array(faces.flatten(), dtype=wp.int32, device=DEV))
    m_i, m_F, m_O, cos_i = t3.mc_flux_3d(mesh, verts, faces, areas, geo, P,
                                         n_ion=12000, n_neu=12000, seed=0)
    floor = (np.abs(cent[:, 0] - 5) < 1.5) & (cent[:, 2] < 9.0)
    nf = int(floor.sum())

    # 2. a TARGET: the floor rate that a recipe of ion energy E*=140 eV would give
    E_true = 140.0
    target = floor_rate_np(m_i, m_F, m_O, cos_i, floor, E_true)
    print(f"floor faces: {nf}")
    print(f"TARGET floor etch rate = {target:.4f}  (the rate a recipe with E*={E_true:.0f} eV gives)")
    print(f"start from E = 100 eV (floor rate {floor_rate_np(m_i,m_F,m_O,cos_i,floor,100.0):.4f})\n")

    # 3. Warp arrays (floor mask as float; only floor faces contribute to the loss)
    def wa(x):
        return wp.array(x.astype(np.float32), dtype=float, device=DEV)
    a_mi, a_mF, a_mO, a_ci = wa(m_i), wa(m_F), wa(m_O), wa(cos_i)
    a_fl = wa(floor.astype(np.float32))
    E = wp.array(np.array([100.0], np.float32), dtype=float, device=DEV, requires_grad=True)

    # 4. gradient descent: optimize E so the floor rate hits the target
    #    (loss is flat in E -> grad ~1e-3, so the step size is large; gradient self-limits at min)
    lr = 18000.0
    for it in range(16):
        sumV = wp.zeros(1, dtype=float, device=DEV, requires_grad=True)
        loss = wp.zeros(1, dtype=float, device=DEV, requires_grad=True)
        tape = wp.Tape()
        with tape:
            wp.launch(floor_rate_sum, dim=len(faces), device=DEV,
                      inputs=[a_mi, a_mF, a_mO, a_ci, a_fl, E, sumV])
            wp.launch(mean_loss, dim=1, device=DEV, inputs=[sumV, float(nf), float(target), loss])
        tape.backward(loss=loss)
        g = float(E.grad.numpy()[0])
        L = float(loss.numpy()[0])
        Eval = float(E.numpy()[0])
        if it % 2 == 0 or it == 15:
            print(f"  iter {it:2d}:  E = {Eval:6.2f} eV   floor_rate = {float(sumV.numpy()[0])/nf:.4f}   "
                  f"loss = {L:.2e}   dL/dE = {g:.2e}")
        # gradient step (autodiff gradient from wp.Tape), clamp to a physical range
        Enew = min(max(Eval - lr * g, 16.0), 400.0)
        E = wp.array(np.array([Enew], np.float32), dtype=float, device=DEV, requires_grad=True)

    print(f"\nRECOVERED recipe: E = {float(E.numpy()[0]):.1f} eV  (true target was {E_true:.0f} eV)")
    print("=> the simulator inverted a target outcome into a recipe, by autodiff. Same machinery")
    print("   scales to many parameters + full-profile targets (inverse design).")


if __name__ == "__main__":
    main()

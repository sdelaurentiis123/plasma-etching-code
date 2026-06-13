"""Neutral transport by radiosity (deterministic, all-bounces) — the ion/neutral split.

Replaces the stochastic many-bounce neutral Monte Carlo with the equilibrium of the diffuse
re-emission integral:

    Gamma = D + (1-s) A Gamma          (A[i,j] = form factor j->i, D = direct sky flux)
    => Gamma = (I - (1-s) A)^-1 D
    adsorbed flux  m_F = s * Gamma      (open-field => Gamma=1 => m_F=s, matches the MC scale)

2D form factor between segments i,j (Lambert, midpoint approximation):
    F_{j->i} = cos_i cos_j / (2 r) * L_i * V_ij
with V_ij the mutual visibility (no other segment blocks the midpoint line) and cos>0 backface
culling. Captures ALL bounces in one linear solve — no 1/sqrt(N) noise, no bounce truncation.

Ions stay on Monte Carlo (directional, single-bounce); only the diffuse neutrals use this.
"""
import numpy as np
from numba import njit


@njit(cache=True, fastmath=True)
def _form_factors(mid, nrm, L, seg):
    """Assemble A[i,j] = F_{j->i} and rowsum_out[i] = sum_j F_{i->j} (for the sky view factor)."""
    M = mid.shape[0]
    A = np.zeros((M, M))
    rowsum_out = np.zeros(M)
    for i in range(M):
        pix = mid[i, 0]; piy = mid[i, 1]; nix = nrm[i, 0]; niy = nrm[i, 1]
        for j in range(M):
            if i == j:
                continue
            dx = mid[j, 0] - pix; dy = mid[j, 1] - piy
            r2 = dx * dx + dy * dy
            if r2 < 1e-12:
                continue
            r = np.sqrt(r2)
            ux = dx / r; uy = dy / r
            cos_i = nix * ux + niy * uy            # angle at i toward j
            cos_j = -(nrm[j, 0] * ux + nrm[j, 1] * uy)  # angle at j toward i
            if cos_i <= 0.0 or cos_j <= 0.0:       # backface cull (must face each other)
                continue
            # visibility: any third segment crossing the i-j midpoint line?
            blocked = False
            for k in range(M):
                if k == i or k == j:
                    continue
                ax = seg[k, 0]; ay = seg[k, 1]; bx = seg[k, 2]; by = seg[k, 3]
                ex = bx - ax; ey = by - ay
                den = ux * ey - uy * ex
                if abs(den) < 1e-14:
                    continue
                t = ((ax - pix) * ey - (ay - piy) * ex) / den
                u = ((ax - pix) * uy - (ay - piy) * ux) / den
                if t > 1e-4 and t < r - 1e-4 and u >= -1e-6 and u <= 1.0 + 1e-6:
                    blocked = True
                    break
            if blocked:
                continue
            ff = cos_i * cos_j / (2.0 * r)
            A[i, j] = ff * L[i]            # F_{j->i}: receiver length L_i
            rowsum_out[i] += ff * L[j]     # F_{i->j}: receiver length L_j
    return A, rowsum_out


def neutral_radiosity(seg, mid, nrm, L, s):
    """Equilibrium adsorbed neutral flux per segment (open-field normalized to s, like MC)."""
    M = len(mid)
    if M == 0:
        return np.zeros(0)
    A, rowsum_out = _form_factors(mid, nrm, L, seg)
    D = np.maximum(0.0, 1.0 - rowsum_out)          # direct sky view factor (reciprocity)
    sys = np.eye(M) - (1.0 - s) * A
    try:
        Gamma = np.linalg.solve(sys, D)
    except np.linalg.LinAlgError:
        Gamma = D
    Gamma = np.clip(Gamma, 0.0, None)
    return s * Gamma                                 # adsorbed flux = m_F

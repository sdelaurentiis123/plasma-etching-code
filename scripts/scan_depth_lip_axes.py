import importlib.util, numpy as np
spec = importlib.util.spec_from_file_location("solve", "scripts/ion_channel_model_solve.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
base = m.Config('barklund_yield', 'b1.7_f0=1', 'zbl', 'unity')
with base:
    ref, _ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
print(f"reference front rate {ref:.3f} nm/s\n")
print(f"{'oxide assignment':<22s} {'energy':<10s} {'C5 depth factor':>16s} {'C6 lip nm/s':>12s}")
for ox in m.OXIDE_OPTIONS:
    for en in m.ENERGY_OPTIONS:
        cfg = m.Config('barklund_yield', ox, en, 'unity')
        with cfg:
            r, _ = m.coupled_rate(1.0, m.ION_SRC, m.E_FRONT)
            lip = m.lip_growth(steps=400)
        print(f"{ox:<22s} {en:<10s} {r/ref:16.3f} {lip:12.3f}")
print("\nC5 required band 0.735 - 0.812 ; C6 required band 0.30 - 0.60 nm/s")

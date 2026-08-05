import importlib.util
spec = importlib.util.spec_from_file_location("solve", "scripts/ion_channel_model_solve.py")
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
print(f"{'oxide assignment':<22s} {'C1a dyn range':>14s} {'C2 peak/normal':>15s}")
for ox in m.OXIDE_OPTIONS:
    cfg = m.Config('barklund_yield', ox, 'zbl', 'unity')
    with cfg:
        dyn, _ = m.gray_curve()
        pk = m.oxide_angular_peak()
    print(f"{ox:<22s} {dyn:14.3f} {pk:15.3f}")
print("\nC1a required 0.20-0.30 (Gray 1993) ; C2 required 1.28-1.36 (Cho/Schaepkens)")
# C6 across the polymer axis, to show the lip is insensitive to the whole space
print()
for po in m.POLYMER_OPTIONS:
    cfg = m.Config(po, 'b1.7_f0=1', 'zbl', 'unity')
    with cfg:
        lip = m.lip_growth(steps=400)
    print(f"C6 lip growth, polymer={po:<22s} {lip:7.3f} nm/s")

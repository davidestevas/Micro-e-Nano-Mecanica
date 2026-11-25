import numpy as np
from pathlib import Path

base = Path(__file__).resolve().parents[1] / "01_lammps_bulk"

results = {}
for label, fname in [("dt = 1.0 fs", "nve_300K_dt1p0.dat"), ("dt = 0.5 fs", "nve_300K_dt0p5.dat")]:
    d = np.loadtxt(base / fname, comments="#")
    t, T, pe, ke, etot = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]
    vcm_max = np.max(np.abs(d[:, 5:8]))

    e0 = etot[0]
    rel_drift_end = abs(etot[-1] - e0) / abs(e0)
    slope = np.polyfit(t, etot, 1)[0]
    rel_drift_fit = abs(slope * (t[-1] - t[0])) / abs(e0)
    std_E = np.std(etot, ddof=1)
    H = np.var(etot, ddof=1) / np.var(ke, ddof=1)

    results[label] = std_E
    print(f"=== NVE 300 K ({label}) ===")
    print(f"  T media                = {np.mean(T):.2f} +/- {np.std(T):.2f} K")
    print(f"  E0                     = {e0:.6f} eV")
    print(f"  std(Etot)              = {std_E:.6e} eV")
    print(f"  Drift rel. (inicio-fim)= {rel_drift_end:.3e}  (criterio < 1e-4: {'OK' if rel_drift_end < 1e-4 else 'FAIL'})")
    print(f"  Drift rel. (ajuste lin)= {rel_drift_fit:.3e}  (criterio < 1e-4: {'OK' if rel_drift_fit < 1e-4 else 'FAIL'})")
    print(f"  Razao H (Eq. 30.8)     = {H:.3e}  (criterio <= 0.05: {'OK' if H <= 0.05 else 'FAIL'})")
    print(f"  max |v_cm|             = {vcm_max:.3e} A/ps")

ratio = results["dt = 1.0 fs"] / results["dt = 0.5 fs"]
p_obs = np.log(ratio) / np.log(2.0)
print(f"\n=== Convergencia Temporal ===")
print(f"  std_E(1.0 fs) / std_E(0.5 fs) = {ratio:.2f}  (esperado ~ 4.00)")
print(f"  Ordem observada p             = {p_obs:.2f}  (formal velocity-Verlet: 2.00)")

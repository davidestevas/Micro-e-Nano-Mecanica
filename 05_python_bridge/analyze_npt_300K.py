import numpy as np
from pathlib import Path

dat = Path(__file__).resolve().parents[1] / "01_lammps_bulk" / "npt_300K_8x8x8.dat"
d = np.loadtxt(dat, comments="#")
t, T, P, a, pe = d[:, 0], d[:, 1], d[:, 2], d[:, 3], d[:, 4]

# Descartar os primeiros 10 ps de transiente
mask = t >= 10.0
t_p, T_p, P_p, a_p = t[mask], T[mask], P[mask], a[mask]

a_mean = np.mean(a_p)
slope_a = np.polyfit(t_p, a_p, 1)[0]
drift_a = slope_a * (t_p[-1] - t_p[0])

print("=== NPT 300 K (8x8x8, producao t = 10..40 ps) ===")
print(f"  <T>          = {np.mean(T_p):.2f} +/- {np.std(T_p, ddof=1):.2f} K")
print(f"  <P>          = {np.mean(P_p):.4f} +/- {np.std(P_p, ddof=1):.4f} GPa")
print(f"  <a(300 K)>   = {a_mean:.6f} A  (expansao vs 0 K: +{(a_mean/3.614925 - 1)*100:.3f}%)")
print(f"  Drift linear = {drift_a:+.3e} A ao longo de 30 ps")

print("\n--- Block averaging para a(300 K) (Eq. 30.20) ---")
dt_samp = t_p[1] - t_p[0]
for nb in [60, 30, 15, 10, 6, 5]:
    lb = len(a_p) // nb
    blocks = a_p[:nb * lb].reshape(nb, lb).mean(axis=1)
    se = np.std(blocks, ddof=1) / np.sqrt(nb)
    print(f"  N_b = {nb:2d} blocos (l_b = {lb * dt_samp:5.2f} ps): SE = {se:.6f} A  (95% CI: +/- {1.96*se:.6f} A)")

import pandas as pd
from pathlib import Path

base = Path(__file__).resolve().parents[1] / "02_lammps_gb_310"
surf = pd.read_csv(base / "surf310_summary.csv")
df = pd.read_csv(base / "gb310_grid_400.csv")

print("=== Superficie livre Cu (310) ===")
for _, r in surf.iterrows():
    print(f"  {r['slab']:>3s}: N={int(r['N_atoms']):4d}, gamma_s = {r['gamma_s_Jm2']:.6f} J/m^2 ({r['gamma_s_Jm2']*1000:.2f} mJ/m^2)")
dg_surf = abs(surf.iloc[1]["gamma_s_Jm2"] - surf.iloc[0]["gamma_s_Jm2"]) * 1000.0
print(f"  Diferenca 40p vs 80p = {dg_surf:.3f} mJ/m^2 (criterio < 1 mJ/m^2: {'OK' if dg_surf < 1.0 else 'CHECK'})")

print(f"\n=== Grelha 10x10x4 Sigma5(310)[001] ({len(df)} configuracoes) ===")
print(f"  max |fz_top| apos minimizacao = {df['fz_top_eVA'].abs().max():.3e} eV/A")
top10 = df.sort_values("gamma_GB_Jm2").head(10)
print("\nTop 10 menores energias:")
print(top10[["ix", "iy", "cf_frac", "dx_A", "dy_A", "N_total", "gamma_GB_Jm2", "fz_top_eVA"]].to_string(index=False))

best = top10.iloc[0]
gs = float(surf.iloc[-1]["gamma_s_Jm2"])
ggb = float(best["gamma_GB_Jm2"])
wsep = 2.0 * gs - ggb
print(f"\n=== Minimo global provisorio (H_grain = 40 A) ===")
print(f"  ix={int(best['ix'])}, iy={int(best['iy'])}, cutoff={best['cf_frac']}*d_nn (rc={best['rc_A']:.4f} A)")
print(f"  dx={best['dx_A']:.5f} A, dy={best['dy_A']:.5f} A, N_total={int(best['N_total'])}")
print(f"  gamma_GB(310) = {ggb:.6f} J/m^2 ({ggb*1000:.2f} mJ/m^2)")
print(f"  gamma_s(310)  = {gs:.6f} J/m^2 ({gs*1000:.2f} mJ/m^2)")
print(f"  Trabalho de separacao W_sep = 2*gamma_s - gamma_GB = {wsep:.6f} J/m^2")

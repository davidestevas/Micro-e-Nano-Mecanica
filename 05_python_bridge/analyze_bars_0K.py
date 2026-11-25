import json
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

trapz = getattr(np, "trapezoid", np.trapz)

ROOT = Path(__file__).resolve().parents[1]
BAR_DIR = ROOT / "04_lammps_bars"
BRIDGE_DIR = ROOT / "05_python_bridge"
FIG_DIR = ROOT / "07_report_figs"

elast = json.loads((BRIDGE_DIR / "cu_elastic_summary.json").read_text())
czm0 = json.loads((BRIDGE_DIR / "czm_0K_params.json").read_text())
M_0K = elast["T_0K"]["M_free_shear_310_GPa"]

summary = {}
fig, ax = plt.subplots(figsize=(8.0, 4.8), dpi=180)

hdr = f"{'L_nom':<7} {'L_mob(A)':>9} {'t_pk(GPa)':>10} {'u_pk(A)':>8} {'maxDrop':>9} {'u_snap(A)':>9} {'dE_end':>8} {'Int_Tn':>8} {'Snap?':>6}"
print(hdr)
print("-" * len(hdr))

for Lnm in [2, 5, 10, 20, 40]:
    p1 = BAR_DIR / f"bar0K_L{Lnm}nm.dat"
    h = p1.read_text().splitlines()[0]
    Lmob0 = float(re.search(r"Lmob0=([0-9.\-]+)", h).group(1))
    Area = float(re.search(r"Area=([0-9.\-]+)", h).group(1))
    d1 = np.loadtxt(p1, comments="#")
    u, tn, dE = d1[:, 1], d1[:, 2], d1[:, 5]

    p2 = BAR_DIR / f"bar0K_L{Lnm}nm_cont.dat"
    if p2.exists():
        d2 = np.loadtxt(p2, comments="#")
        u = np.concatenate([u, d2[:, 1]])
        tn = np.concatenate([tn, d2[:, 2]])
        dE = np.concatenate([dE, d2[:, 5]])

    # Guardar curva combinada limpa para comparacao direta com o Abaqus
    np.savetxt(BAR_DIR / f"bar0K_L{Lnm}nm_full.dat", np.column_stack([u, tn, dE]),
               header=f"L_nom={Lnm} nm | Lmob0={Lmob0:.6f} A | Area={Area:.6f} A^2\n# u_grip_A tn_avg_GPa dE_Jm2")

    idx_pk = int(np.argmax(tn))
    t_pk = float(tn[idx_pk])
    u_pk = float(u[idx_pk])
    drops = -np.diff(tn)
    idx_drop = int(np.argmax(drops))
    max_drop = float(drops[idx_drop])
    u_snap = float(u[idx_drop + 1])
    W_int = float(trapz(tn, u) * 0.1)
    snap = max_drop > 5.0

    print(f"{Lnm:2d} nm   {Lmob0:9.3f} {t_pk:10.4f} {u_pk:8.2f} {max_drop:9.3f} {u_snap:9.2f} {dE[-1]:8.4f} {W_int:8.4f} {'SIM' if snap else 'NAO':>6}")
    summary[f"L{Lnm}nm"] = {
        "L_nom_nm": Lnm, "L_mob_A": Lmob0, "Area_A2": Area,
        "t_peak_GPa": t_pk, "u_peak_A": u_pk, "u_snap_A": u_snap,
        "max_drop_GPa": max_drop, "dE_end_Jm2": float(dE[-1]),
        "W_int_Jm2": W_int, "snap_back": snap,
    }

    # Converter para forca em nN (1 GPa * 1 A^2 = 0.01 nN) e abertura em nm
    F_nN = tn * Area * 0.01
    ax.plot(u * 0.1, tn, lw=1.7, label=rf"$L = {Lnm}$ nm ($t_{{\max}}={t_pk:.2f}$ GPa, $u_{{\mathrm{{pk}}}}={u_pk*0.1:.2f}$ nm)")

(BRIDGE_DIR / "bars_0K_summary.json").write_text(json.dumps(summary, indent=2))

ax.axhline(0, color="k", lw=0.6, ls=":")
ax.set_xlabel("Deslocamento total imposto $u$ (nm)")
ax.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax.set_title(r"Barras bicristal Cu $\Sigma 5(310)$ a 0 K (LAMMPS MS): efeito do comprimento $L$")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8.5)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_bars_0K_lammps.png")
print(f"\nResumo guardado em: {BRIDGE_DIR / 'bars_0K_summary.json'}")
print(f"Figura guardada em: {FIG_DIR / 'fig_bars_0K_lammps.png'}")

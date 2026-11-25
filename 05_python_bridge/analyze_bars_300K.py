import json
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
BAR_DIR = ROOT / "04_lammps_bars"
BRIDGE_DIR = ROOT / "05_python_bridge"
FIG_DIR = ROOT / "07_report_figs"

def moving_avg(y, win=15):
    if len(y) < win:
        return y.copy()
    pad = win // 2
    yp = np.pad(y, (pad, pad), mode="edge")
    return np.convolve(yp, np.ones(win) / win, mode="valid")[:len(y)]

def last_frame_gap(dump_path):
    lines = dump_path.read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    blk = lines[starts[-1]:]
    ai = [i for i, l in enumerate(blk) if l.startswith("ITEM: ATOMS")][0]
    pts = np.array([[float(x) for x in r.split()[:5]] for r in blk[ai + 1:] if r.strip()])
    zs = np.sort(pts[:, 4])
    gaps = np.diff(zs)
    ig = int(np.argmax(gaps))
    return float(gaps[ig]), float(zs[ig]), float(zs[ig + 1])

tags = ["L5nm_v10ms", "L5nm_v10ms_s2", "L5nm_v1ms", "L10nm_v10ms", "L20nm_v10ms", "L40nm_v10ms"]
summary = {}

hdr = f"{'Tag':<15} {'Lmob(A)':>8} {'t_pk_sm':>8} {'t_pk_rw':>8} {'u_pk(A)':>8} {'u_end':>6} {'tn_end':>8} {'gap_end':>8} {'Estado':>10}"
print(hdr)
print("-" * len(hdr))

fig, ax = plt.subplots(figsize=(8.6, 5.0), dpi=180)

for tag in tags:
    meta = (BAR_DIR / f"bar300K_{tag}_meta.txt").read_text()
    Lmob0 = float(re.search(r"Lmob0=([0-9.\-]+)", meta).group(1))
    Area = float(re.search(r"Area=([0-9.\-]+)", meta).group(1))
    vel = float(re.search(r"vel=([0-9.\-]+)", meta).group(1))

    d1 = np.loadtxt(BAR_DIR / f"bar300K_{tag}.dat", comments="#")
    u, tn = d1[:, 1], d1[:, 2]
    dmp_p = BAR_DIR / f"dump.bar300K_{tag}.lammpstrj"

    cont_p = BAR_DIR / f"bar300K_{tag}_cont.dat"
    if cont_p.exists():
        d2 = np.loadtxt(cont_p, comments="#")
        u0_c = float(d2[0, 1])
        keep = u < u0_c
        u = np.concatenate([u[keep], d2[:, 1]])
        tn = np.concatenate([tn[keep], d2[:, 2]])
        dmp_p = BAR_DIR / f"dump.bar300K_{tag}_cont.lammpstrj"

    np.savetxt(BAR_DIR / f"bar300K_{tag}_full.dat", np.column_stack([u, tn]),
               header=f"tag={tag} | Lmob0={Lmob0:.6f} A | Area={Area:.6f} A^2 | vel={vel} A/ps\n# u_grip_A tn_avg_GPa")

    tn_s = moving_avg(tn, 15)
    idx_pk = int(np.argmax(tn_s))
    t_pk_sm = float(tn_s[idx_pk])
    t_pk_rw = float(np.max(tn))
    u_pk = float(u[idx_pk])
    tn_end_avg = float(np.mean(tn[-30:]))

    gap, zg1, zg2 = last_frame_gap(dmp_p)
    broke = gap > 3.5

    print(f"{tag:<15} {Lmob0:8.3f} {t_pk_sm:8.3f} {t_pk_rw:8.3f} {u_pk:8.2f} {u[-1]:6.2f} {tn_end_avg:+8.3f} {gap:7.2f}A {'PARTIU' if broke else 'INTACTA':>10}")

    summary[tag] = {
        "Lmob0_A": Lmob0, "Area_A2": Area, "vel_A_per_ps": vel,
        "t_peak_smooth_GPa": t_pk_sm, "t_peak_raw_GPa": t_pk_rw,
        "u_peak_A": u_pk, "u_end_A": float(u[-1]),
        "tn_end_avg_GPa": tn_end_avg, "max_gap_end_A": gap,
        "gap_z_range_A": [zg1, zg2], "broke": broke,
    }

    ls = "--" if ("s2" in tag or "v1ms" in tag) else "-"
    ax.plot(u * 0.1, tn_s, ls, lw=1.6, label=rf"{tag} ($t_{{\max}}={t_pk_sm:.2f}$ GPa, $u_{{\mathrm{{pk}}}}={u_pk*0.1:.2f}$ nm)")

(BRIDGE_DIR / "bars_300K_summary.json").write_text(json.dumps(summary, indent=2))

ax.axhline(0, color="k", lw=0.6, ls=":")
ax.set_xlabel("Deslocamento total imposto $u$ (nm)")
ax.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax.set_title(r"Barras bicristal Cu $\Sigma 5(310)$ a 300 K (LAMMPS MD)")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_bars_300K_lammps.png")
print(f"\nFigura guardada em: {FIG_DIR / 'fig_bars_300K_lammps.png'}")

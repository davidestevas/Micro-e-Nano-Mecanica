import json
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

trapz = getattr(np, "trapezoid", np.trapz)

ROOT = Path(__file__).resolve().parents[1]
PULL_DIR = ROOT / "03_lammps_czm_pull"
BRIDGE_DIR = ROOT / "05_python_bridge"
FIG_DIR = ROOT / "07_report_figs"

elast = json.loads((BRIDGE_DIR / "cu_elastic_summary.json").read_text())
czm0 = json.loads((BRIDGE_DIR / "czm_0K_params.json").read_text())

M_0K = elast["T_0K"]["M_free_shear_310_GPa"]
M_300K = elast["T_300K"]["M_free_shear_310_GPa"]
E_300K = elast["T_300K"]["E_310_GPa"]

def moving_avg(y, win=15):
    pad = win // 2
    yp = np.pad(y, (pad, pad), mode="edge")
    return np.convolve(yp, np.ones(win) / win, mode="valid")[:len(y)]

def xu_needleman(delta, t0, dc):
    x = np.maximum(delta, 0.0) / dc
    return t0 * x * np.exp(1.0 - x)

def load_300K(vtag):
    meta = (PULL_DIR / f"pull300K_{vtag}_meta.txt").read_text()
    Lmob0 = float(re.search(r"Lmob0=([0-9.\-]+)", meta).group(1))
    Area = float(re.search(r"Area=([0-9.\-]+)", meta).group(1))
    d = np.loadtxt(PULL_DIR / f"pull300K_{vtag}.dat", comments="#")
    # step, u_grip, tn_avg, tn_top, tn_bot, ts_x, ts_y, T_xy, T_all
    return {
        "vtag": vtag, "Lmob0": Lmob0, "Area": Area,
        "u": d[:, 1], "tn": d[:, 2], "tn_top": d[:, 3], "tn_bot": d[:, 4],
        "tsx": d[:, 5], "tsy": d[:, 6], "T_xy": d[:, 7], "T_all": d[:, 8],
    }

# Carregar referencia 0 K (Hmob11)
d0_raw = np.loadtxt(PULL_DIR / "pull0K_Hmob11.dat", comments="#")
h0 = (PULL_DIR / "pull0K_Hmob11.dat").read_text().splitlines()[0]
L0_ref = float(re.search(r"Lmob0=([0-9.\-]+)", h0).group(1))
d_0K = d0_raw[:, 1] - d0_raw[:, 2] * L0_ref / M_0K
tn_0K = d0_raw[:, 2]

results = {}
for vtag, win in [("v10ms", 15), ("v1ms", 15)]:
    p = PULL_DIR / f"pull300K_{vtag}.dat"
    if not p.exists():
        continue
    c = load_300K(vtag)
    tn_s = moving_avg(c["tn"], win)
    delta = c["u"] - c["tn"] * c["Lmob0"] / M_300K
    delta_s = c["u"] - tn_s * c["Lmob0"] / M_300K

    idx_pk = int(np.argmax(tn_s))
    t0_s = float(tn_s[idx_pk])
    t0_raw = float(np.max(c["tn"]))
    d_pk = float(delta_s[idx_pk])

    # Integrar trabalho de separacao (ate separacao completa, sem ringing final)
    W_tot = float(trapz(c["tn"], delta) * 0.1)
    W_u = float(trapz(c["tn"], c["u"]) * 0.1)
    rms_imb = float(np.sqrt(np.mean((c["tn_top"] - c["tn_bot"]) ** 2)))

    # Ajuste bilinear (conserva t0_s e W_tot)
    mask_el = (np.arange(len(tn_s)) < idx_pk) & (tn_s >= 0.10 * t0_s) & (tn_s <= 0.55 * t0_s)
    K0 = float(np.polyfit(delta_s[mask_el], tn_s[mask_el], 1)[0])
    d0_bil = t0_s / K0
    df_bil = 2.0 * (W_tot / 0.1) / t0_s

    # Ajuste Xu-Needleman (exato em t0_s, W_tot e LSQ livre)
    dc_exact = (W_tot / 0.1) / (np.e * t0_s)
    popt, _ = curve_fit(xu_needleman, np.maximum(delta_s, 0.0), tn_s, p0=[t0_s, dc_exact], bounds=([1.0, 0.05], [40.0, 5.0]))
    t0_xn, dc_xn = float(popt[0]), float(popt[1])
    Gc_xn = float(np.e * t0_xn * dc_xn * 0.1)

    l_cz_M = (M_300K * 1e9) * W_tot / (t0_s * 1e9) ** 2 * 1e9
    l_cz_E = (E_300K * 1e9) * W_tot / (t0_s * 1e9) ** 2 * 1e9
    L_snap_bil_M = M_300K * (df_bil - d0_bil) / t0_s * 0.1
    L_snap_bil_E = E_300K * (df_bil - d0_bil) / t0_s * 0.1

    print(f"\n=== 300 K ({vtag}, L_mob = {c['Lmob0']:.3f} A) ===")
    print(f"  <T_xy> = {np.mean(c['T_xy']):.1f} +/- {np.std(c['T_xy']):.1f} K | <T_all> = {np.mean(c['T_all']):.1f} +/- {np.std(c['T_all']):.1f} K")
    print(f"  RMS |top - bot|        = {rms_imb:.3f} GPa | tn_final = {np.mean(c['tn'][-50:]):+.4f} GPa")
    print(f"  Pico tn^0 (suav / raw) = {t0_s:.3f} / {t0_raw:.3f} GPa (em delta = {d_pk:.3f} A)")
    print(f"  Energia G_c (d / u)    = {W_tot:.4f} / {W_u:.4f} J/m^2")
    print(f"  Bilinear               : t0={t0_s:.3f} GPa, Gc={W_tot:.4f} J/m^2, K0={K0:.2f} GPa/A, d0={d0_bil:.3f} A, df={df_bil:.3f} A")
    print(f"  Xu-Needleman (exata)   : t0={t0_s:.3f} GPa, dc={dc_exact:.4f} A, Gc={W_tot:.4f} J/m^2")
    print(f"  Xu-Needleman (LSQ)     : t0={t0_xn:.3f} GPa, dc={dc_xn:.4f} A, Gc={Gc_xn:.4f} J/m^2")
    print(f"  l_cz (M / E)           = {l_cz_M:.3f} / {l_cz_E:.3f} nm | L_snap Bilinear (PBC / livre) = {L_snap_bil_M:.3f} / {L_snap_bil_E:.3f} nm")

    results[vtag] = {
        "Lmob0_A": c["Lmob0"], "t0_smooth_GPa": t0_s, "t0_raw_GPa": t0_raw,
        "Gc_Jm2": W_tot, "K0_GPa_per_A": K0, "d0_bil_A": d0_bil, "df_bil_A": df_bil,
        "dc_xn_exact_A": dc_exact, "t0_xn_lsq_GPa": t0_xn, "dc_xn_lsq_A": dc_xn, "Gc_xn_lsq_Jm2": Gc_xn,
        "l_cz_M_nm": l_cz_M, "l_cz_E_nm": l_cz_E,
        "L_snap_bil_PBC_nm": L_snap_bil_M, "L_snap_bil_free_nm": L_snap_bil_E,
        "_delta_s": delta_s, "_tn_s": tn_s, "_tn_raw": c["tn"],
    }

# Guardar JSON limpo
export = {k: {kk: vv for kk, vv in v.items() if not kk.startswith("_")} for k, v in results.items()}
(BRIDGE_DIR / "czm_300K_params.json").write_text(json.dumps(export, indent=2))

# Figura comparativa 0 K vs 300 K (10 m/s e 1 m/s)
fig, ax = plt.subplots(figsize=(7.2, 4.6), dpi=180)
ax.plot(d_0K, tn_0K, "k-", lw=2.0, label=rf"0 K (MS, $t_0={czm0['bilinear_0K']['t0_GPa']:.2f}$ GPa, $G_c={czm0['bilinear_0K']['Gc_Jm2']:.2f}$ J/m$^2$)")
for vtag, col, lbl in [("v10ms", "C0", "300 K (10 m/s)"), ("v1ms", "C1", "300 K (1 m/s)")]:
    if vtag in results:
        r = results[vtag]
        ax.plot(r["_delta_s"], r["_tn_s"], color=col, lw=1.8,
                label=rf"{lbl}: $t_0={r['t0_smooth_GPa']:.2f}$ GPa, $G_c={r['Gc_Jm2']:.2f}$ J/m$^2$")
if "v1ms" in results:
    r = results["v1ms"]
    dd = np.linspace(0, 6.0, 400)
    ax.plot([0, r["d0_bil_A"], r["df_bil_A"]], [0, r["t0_smooth_GPa"], 0], "C1--", lw=1.6, label="Bilinear (300 K, 1 m/s)")
    ax.plot(dd, xu_needleman(dd, r["t0_smooth_GPa"], r["dc_xn_exact_A"]), "C3:", lw=2.0, label="Xu-Needleman (300 K, 1 m/s)")

ax.axhline(0, color="k", lw=0.6, ls=":")
ax.set_xlim(-0.2, 5.5)
ax.set_xlabel(r"Separação coesiva $\delta$ ($\mathrm{\AA}$)")
ax.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax.set_title(r"Lei coesiva Cu $\Sigma 5(310)[001]$: efeito da temperatura e taxa")
ax.grid(True, alpha=0.3)
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(FIG_DIR / "fig_czm_0K_vs_300K.png")
print(f"\nFigura guardada em: {FIG_DIR / 'fig_czm_0K_vs_300K.png'}")

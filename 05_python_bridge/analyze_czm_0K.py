import json
import re
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

trapz = getattr(np, "trapezoid", np.trapz)

ROOT = Path(__file__).resolve().parents[1]
PULL_DIR = ROOT / "03_lammps_czm_pull"
FIG_DIR = ROOT / "07_report_figs"
BRIDGE_DIR = ROOT / "05_python_bridge"

elast = json.loads((BRIDGE_DIR / "cu_elastic_summary.json").read_text())
C33 = elast["T_0K"]["C3333_310_GPa"]          # 188.84 GPa
Mfs = elast["T_0K"]["M_free_shear_310_GPa"]   # 177.73 GPa
E310 = elast["T_0K"]["E_310_GPa"]             # 81.47 GPa

W_sep_ref = 2.0 * 1.515643 - 0.904776         # 2.126510 J/m^2

def load_pull(tag):
    p = PULL_DIR / f"pull0K_{tag}.dat"
    header = p.read_text().splitlines()[0]
    Lmob0 = float(re.search(r"Lmob0=([0-9.\-]+)", header).group(1))
    Area = float(re.search(r"Area=([0-9.\-]+)", header).group(1))
    d = np.loadtxt(p, comments="#")
    return {
        "tag": tag, "Lmob0": Lmob0, "Area": Area,
        "u": d[:, 1], "tn": d[:, 2], "tn_top": d[:, 3], "tn_bot": d[:, 4],
        "tsx": d[:, 5], "tsy": d[:, 6], "dE": d[:, 7],
    }

def xu_needleman(delta, t0, dc):
    x = np.maximum(delta, 0.0) / dc
    return t0 * x * np.exp(1.0 - x)

tags = [t for t in ["Hmob05", "Hmob08", "Hmob11", "Hmob15", "Hmob25"] if (PULL_DIR / f"pull0K_{t}.dat").exists()]
cases = {}

print(f"Referencia Passo 2: W_sep = 2*gamma_s(310) - gamma_GB(310) = {W_sep_ref:.6f} J/m^2\n")
hdr = f"{'Caso':<8} {'L_mob(A)':>8} {'t0(GPa)':>8} {'u_pk(A)':>7} {'d_pk(A)':>7} {'maxDrop':>8} {'dE/A':>8} {'Int_Tn':>8} {'Err_Int':>8} {'Check2%':>8}"
print(hdr)
print("-" * len(hdr))

passing = []
for tag in tags:
    c = load_pull(tag)
    d_mfs = c["u"] - c["tn"] * c["Lmob0"] / Mfs
    c["d_mfs"] = d_mfs

    idx_pk = int(np.argmax(c["tn"]))
    t0 = float(c["tn"][idx_pk])
    u_pk = float(c["u"][idx_pk])
    d_pk = float(d_mfs[idx_pk])
    max_drop = float(np.max(-np.diff(c["tn"])))
    W_dE = float(c["dE"][-1])
    W_int = float(trapz(c["tn"], d_mfs) * 0.1)
    err_int = 100.0 * (W_int - W_sep_ref) / W_sep_ref
    ok = abs(err_int) <= 2.0

    c.update({"t0": t0, "d_pk": d_pk, "W_dE": W_dE, "W_int": W_int, "err_int": err_int, "max_drop": max_drop, "ok": ok})
    cases[tag] = c
    if ok:
        passing.append(tag)

    print(f"{tag:<8} {c['Lmob0']:8.3f} {t0:8.3f} {u_pk:7.2f} {d_pk:7.4f} {max_drop:8.3f} {W_dE:8.4f} {W_int:8.4f} {err_int:+7.2f}% {'OK' if ok else 'SNAP':>8}")

# Escolher o caso mais espesso que passa o check de 2% (ou o de menor erro)
ref_tag = passing[-1] if passing else min(tags, key=lambda t: abs(cases[t]["err_int"]))
ref = cases[ref_tag]
d = ref["d_mfs"]
tn = ref["tn"]
idx_pk = int(np.argmax(tn))
t0 = ref["t0"]
Gc = ref["W_int"]

# 1. Ajuste Bilinear (conserva t0 e Gc; rigidez inicial por fit linear entre 5% e 50% do pico)
mask_el = (np.arange(len(tn)) < idx_pk) & (tn >= 0.05 * t0) & (tn <= 0.50 * t0)
K0_GPa_A = float(np.polyfit(d[mask_el], tn[mask_el], 1)[0]) if np.sum(mask_el) >= 3 else float(t0 / d[idx_pk])
d0_bil = t0 / K0_GPa_A
df_bil = 2.0 * (Gc / 0.1) / t0

# 2. Ajuste Xu-Needleman: (a) conservando t0 e Gc; (b) minimos quadrados livres
dc_xn_exact = (Gc / 0.1) / (np.e * t0)
popt, _ = curve_fit(xu_needleman, np.maximum(d, 0.0), tn, p0=[t0, dc_xn_exact], bounds=([1.0, 0.05], [50.0, 5.0]))
t0_xn_fit, dc_xn_fit = float(popt[0]), float(popt[1])
Gc_xn_fit = float(np.e * t0_xn_fit * dc_xn_fit * 0.1)

# Declive maximo negativo da curva atomistica (para prever L_snap exato)
post_mask = np.arange(len(d)) >= idx_pk
dd_post = np.diff(d[post_mask])
dtn_post = np.diff(tn[post_mask])
valid = dd_post > 1e-4
max_soft_slope = float(np.max(-dtn_post[valid] / dd_post[valid])) if np.any(valid) else np.nan
L_snap_atom_Mfs = (Mfs / max_soft_slope) * 0.1 if max_soft_slope > 0 else np.nan
L_snap_atom_E = (E310 / max_soft_slope) * 0.1 if max_soft_slope > 0 else np.nan

l_cz_Mfs = (Mfs * 1e9) * Gc / (t0 * 1e9) ** 2 * 1e9
l_cz_E = (E310 * 1e9) * Gc / (t0 * 1e9) ** 2 * 1e9
L_snap_bil_Mfs = Mfs * (df_bil - d0_bil) / t0 * 0.1
L_snap_bil_E = E310 * (df_bil - d0_bil) / t0 * 0.1
L_snap_xn_Mfs = l_cz_Mfs
L_snap_xn_E = l_cz_E

print(f"\n============================================================")
print(f"LEI COESIVA DE REFERENCIA A 0 K: {ref_tag} (L_mob = {ref['Lmob0']:.3f} A)")
print(f"============================================================")
print(f"Erro da area vs W_sep = {ref['err_int']:+.3f}%  (dE/A = {ref['W_dE']:.6f} J/m^2)")
print(f"Bilinear (t0 e Gc exatos):")
print(f"  t0 = {t0:.4f} GPa | Gc = {Gc:.6f} J/m^2 | K0 = {K0_GPa_A:.3f} GPa/A ({K0_GPa_A*1e19:.3e} N/m^3)")
print(f"  d0 = {d0_bil:.4f} A | df = {df_bil:.4f} A")
print(f"Xu-Needleman (t0 e Gc exatos):")
print(f"  t0 = {t0:.4f} GPa | dc = {dc_xn_exact:.4f} A | Gc = {Gc:.6f} J/m^2")
print(f"Xu-Needleman (LSQ livre):")
print(f"  t0 = {t0_xn_fit:.4f} GPa | dc = {dc_xn_fit:.4f} A | Gc = {Gc_xn_fit:.6f} J/m^2 (erro vs W_sep: {(Gc_xn_fit/W_sep_ref-1)*100:+.2f}%)")
print(f"\nDeclive max de amolecimento atomistico |dT/dd|_max = {max_soft_slope:.2f} GPa/A")
print(f"Previsoes de L_snap (0 K):")
print(f"  Curva atomistica : {L_snap_atom_Mfs:.3f} nm (PBC laterais) | {L_snap_atom_E:.3f} nm (barra livre)")
print(f"  Xu-Needleman     : {L_snap_xn_Mfs:.3f} nm (PBC laterais) | {L_snap_xn_E:.3f} nm (barra livre)")
print(f"  Bilinear         : {L_snap_bil_Mfs:.3f} nm (PBC laterais) | {L_snap_bil_E:.3f} nm (barra livre)")

out = {
    "ref_tag": ref_tag,
    "W_sep_ref_Jm2": W_sep_ref,
    "M_free_shear_0K_GPa": Mfs,
    "C3333_0K_GPa": C33,
    "E310_0K_GPa": E310,
    "max_softening_slope_GPa_per_A": max_soft_slope,
    "bilinear_0K": {
        "t0_GPa": t0, "Gc_Jm2": Gc, "K0_GPa_per_A": K0_GPa_A,
        "K0_N_per_m3": K0_GPa_A * 1e19, "d0_A": d0_bil, "df_A": df_bil,
        "L_snap_PBC_nm": L_snap_bil_Mfs, "L_snap_free_nm": L_snap_bil_E,
    },
    "xu_needleman_0K_exact": {
        "t0_GPa": t0, "dc_A": dc_xn_exact, "Gc_Jm2": Gc,
        "L_snap_PBC_nm": L_snap_xn_Mfs, "L_snap_free_nm": L_snap_xn_E,
    },
    "xu_needleman_0K_lsq": {
        "t0_GPa": t0_xn_fit, "dc_A": dc_xn_fit, "Gc_Jm2": Gc_xn_fit,
    },
}
(BRIDGE_DIR / "czm_0K_params.json").write_text(json.dumps(out, indent=2))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.8, 4.5), dpi=180)
for tag in tags:
    c = cases[tag]
    ax1.plot(c["u"], c["tn"], lw=1.6, label=rf"{tag} ($L_{{\mathrm{{mob}}}}={c['Lmob0']:.1f}\ \mathrm{{\AA}}$)")
    ax2.plot(c["d_mfs"], c["tn"], lw=1.5, alpha=0.85, label=rf"{tag} (err {c['err_int']:+.1f}%)")

dd = np.linspace(0, 6.0, 400)
ax2.plot([0, d0_bil, df_bil], [0, t0, 0], "k--", lw=1.8, label=rf"Bilinear ($t_0={t0:.2f}$ GPa)")
ax2.plot(dd, xu_needleman(dd, t0, dc_xn_exact), "C3:", lw=2.2, label="Xu-Needleman (exata)")
ax2.plot(dd, xu_needleman(dd, t0_xn_fit, dc_xn_fit), "C2-.", lw=1.6, label="Xu-Needleman (LSQ)")

ax1.set_xlabel(r"Abertura total $u_{\mathrm{grip}}$ ($\mathrm{\AA}$)")
ax1.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax1.set_title("Curvas brutas a 0 K: transição para snap-back")
ax1.grid(True, alpha=0.3)
ax1.legend(fontsize=8)

ax2.set_xlabel(r"Separação coesiva $\delta = u - T_n L_{\mathrm{mob}}/M$ ($\mathrm{\AA}$)")
ax2.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax2.set_title(rf"Lei coesiva Cu $\Sigma 5(310)$ a 0 K ($W_{{\mathrm{{sep}}}}={W_sep_ref:.3f}$ J/m$^2$)")
ax2.set_xlim(-0.2, 5.5)
ax2.grid(True, alpha=0.3)
ax2.legend(fontsize=7.5)

fig.tight_layout()
fig.savefig(FIG_DIR / "fig_czm_0K_check.png")
print(f"\nFigura guardada em: {FIG_DIR / 'fig_czm_0K_check.png'}")

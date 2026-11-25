import json
import re
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
CZM_DIR = ROOT / "03_lammps_czm_pull"
BAR_DIR = ROOT / "04_lammps_bars"
BRIDGE_DIR = ROOT / "05_python_bridge"
ABQ_DIR = ROOT / "06_abaqus_bars"
FIG_DIR = ROOT / "07_report_figs"
FIG_DIR.mkdir(exist_ok=True)

if not (ABQ_DIR / "bar0K_L5nm_xn_curve.csv").exists():
    sys.exit("ERRO: Ficheiros *_curve.csv do Abaqus nao encontrados em 06_abaqus_bars/. Faz 'git pull' (ou descompacta abaqus_results_for_mac.zip) primeiro.")

plt.rcParams.update({
    "font.size": 9.5,
    "axes.labelsize": 10,
    "axes.titlesize": 10.5,
    "legend.fontsize": 8.0,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "lines.linewidth": 1.6,
    "axes.grid": True,
    "grid.alpha": 0.28,
})

elast = json.loads((BRIDGE_DIR / "cu_elastic_summary.json").read_text())
czm0 = json.loads((BRIDGE_DIR / "czm_0K_params.json").read_text())
czm300 = json.loads((BRIDGE_DIR / "czm_300K_params.json").read_text())
bars0 = json.loads((BRIDGE_DIR / "bars_0K_summary.json").read_text())
bars300 = json.loads((BRIDGE_DIR / "bars_300K_summary.json").read_text())

M0 = elast["T_0K"]["M_free_shear_310_GPa"]
M300 = elast["T_300K"]["M_free_shear_310_GPa"]

def moving_avg(y, win=15):
    if len(y) < win:
        return y.copy()
    pad = win // 2
    yp = np.pad(y, (pad, pad), mode="edge")
    return np.convolve(yp, np.ones(win) / win, mode="valid")[:len(y)]

def load_abq_csv(jobname):
    p = ABQ_DIR / f"{jobname}_curve.csv"
    if not p.exists():
        return None
    d = np.loadtxt(p, delimiter=",", skiprows=1)
    if d.ndim == 1:
        d = d[None, :]
    return {"u": d[:, 3], "delta": d[:, 4], "F": d[:, 5], "Tn": d[:, 6], "ALLDMD": d[:, 9]}

def xu_needleman(delta, t0, dc):
    x = np.maximum(delta, 0.0) / dc
    return t0 * x * np.exp(1.0 - x)

def save_fig(fig, stem):
    fig.savefig(FIG_DIR / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG_DIR / f"{stem}.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

# FIG 1: Leis Coesivas (0 K e 300 K)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
for tag, col, lbl in [
    ("Hmob05", "C4", r"$L_{\mathrm{mob}}=0.45$ nm"),
    ("Hmob08", "C0", r"$L_{\mathrm{mob}}=0.77$ nm (OK)"),
    ("Hmob11", "k",  r"$L_{\mathrm{mob}}=1.06$ nm (Ref.)"),
    ("Hmob15", "C3", r"$L_{\mathrm{mob}}=1.57$ nm (snap)"),
    ("Hmob25", "C1", r"$L_{\mathrm{mob}}=2.55$ nm (snap)"),
]:
    p = CZM_DIR / f"pull0K_{tag}.dat"
    h = p.read_text().splitlines()[0]
    Lm = float(re.search(r"Lmob0=([0-9.\-]+)", h).group(1))
    arr = np.loadtxt(p, comments="#")
    d_nm = (arr[:, 1] - arr[:, 2] * Lm / M0) * 0.1
    ax1.plot(d_nm, arr[:, 2], color=col, lw=2.0 if tag == "Hmob11" else 1.3, label=lbl)

dd = np.linspace(0, 0.55, 400)
b0 = czm0["bilinear_0K"]
x0 = czm0["xu_needleman_0K_exact"]
ax1.plot([0, b0["d0_A"]*0.1, b0["df_A"]*0.1], [0, b0["t0_GPa"], 0], "k--", lw=1.5, label="Bilinear (0 K)")
ax1.plot(dd, xu_needleman(dd, x0["t0_GPa"], x0["dc_A"]*0.1), "C2:", lw=2.1, label="Xu-Needleman (0 K)")
ax1.set_xlim(-0.01, 0.52)
ax1.set_ylim(-0.5, 17.5)
ax1.set_xlabel(r"Separação coesiva $\delta$ (nm)")
ax1.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax1.set_title(r"(a) Extração a 0 K: transição para snap-back ($L_{\mathrm{snap}}=1.12$ nm)")
ax1.legend(loc="upper right", fontsize=7.5)

p0 = CZM_DIR / "pull0K_Hmob11.dat"
Lm0 = float(re.search(r"Lmob0=([0-9.\-]+)", p0.read_text().splitlines()[0]).group(1))
arr0 = np.loadtxt(p0, comments="#")
ax2.plot((arr0[:, 1] - arr0[:, 2]*Lm0/M0)*0.1, arr0[:, 2], "k-", lw=2.0,
         label=rf"0 K MS ($t_0={b0['t0_GPa']:.2f}$ GPa, $G_c={b0['Gc_Jm2']:.2f}$ J/m$^2$)")

for vtag, col, lbl in [("v10ms", "C0", "300 K (10 m/s)"), ("v1ms", "C1", "300 K (1 m/s)")]:
    meta = (CZM_DIR / f"pull300K_{vtag}_meta.txt").read_text()
    Lm = float(re.search(r"Lmob0=([0-9.\-]+)", meta).group(1))
    arr = np.loadtxt(CZM_DIR / f"pull300K_{vtag}.dat", comments="#")
    tn_s = moving_avg(arr[:, 2], 15)
    d_nm = (arr[:, 1] - tn_s * Lm / M300) * 0.1
    r = czm300[vtag]
    ax2.plot(d_nm, tn_s, color=col, lw=1.7,
             label=rf"{lbl} ($t_0={r['t0_smooth_GPa']:.2f}$ GPa, $G_c={r['Gc_Jm2']:.2f}$ J/m$^2$)")

r1 = czm300["v1ms"]
dd2 = np.linspace(0, 0.70, 400)
ax2.plot([0, r1["d0_bil_A"]*0.1, r1["df_bil_A"]*0.1], [0, r1["t0_smooth_GPa"], 0], "C1--", lw=1.4, label="Bilinear (300 K, 1 m/s)")
ax2.plot(dd2, xu_needleman(dd2, r1["t0_smooth_GPa"], r1["dc_xn_exact_A"]*0.1), "C3:", lw=2.0, label="Xu-Needleman (300 K, 1 m/s)")
ax2.set_xlim(-0.01, 0.65)
ax2.set_ylim(-0.5, 17.5)
ax2.set_xlabel(r"Separação coesiva $\delta$ (nm)")
ax2.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax2.set_title("(b) Efeito da temperatura e da taxa de carregamento")
ax2.legend(loc="upper right", fontsize=7.5)
fig.tight_layout()
save_fig(fig, "fig01_czm_laws")

# FIG 2: Transicao para snap-back no Abaqus
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
cols_L = {1: "C0", 2: "C2", 3: "C1", 5: "C3"}
for L in [1, 2, 3, 5]:
    cx = load_abq_csv(f"bar0K_L{L}nm_xn")
    cb = load_abq_csv(f"bar0K_L{L}nm_bilin")
    if cx:
        ax1.plot(cx["u"], cx["Tn"], "-", color=cols_L[L], lw=1.8, label=rf"$L={L}$ nm (XN)")
    if cb:
        ax1.plot(cb["u"], cb["Tn"], "--", color=cols_L[L], lw=1.3, alpha=0.85, label=rf"$L={L}$ nm (Bilin)")
ax1.set_xlim(0, 0.55)
ax1.set_ylim(0, 17.0)
ax1.set_xlabel("Deslocamento total $u$ (nm)")
ax1.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax1.set_title(r"(a) Abaqus Riks a 0 K: $L_{\mathrm{snap}}^{\mathrm{XN}}=1.53$ nm vs. $L_{\mathrm{snap}}^{\mathrm{Bilin}}=2.74$ nm")
ax1.legend(ncol=2, fontsize=7.3)

for L in [2, 3, 5]:
    cx = load_abq_csv(f"bar300K_L{L}nm_xn")
    cb = load_abq_csv(f"bar300K_L{L}nm_bilin")
    if cx:
        ax2.plot(cx["u"], cx["Tn"], "-", color=cols_L[L], lw=1.8, label=rf"$L={L}$ nm (XN, 10 m/s)")
    if cb:
        ax2.plot(cb["u"], cb["Tn"], "--", color=cols_L[L], lw=1.3, alpha=0.85, label=rf"$L={L}$ nm (Bilin, 10 m/s)")
cx1 = load_abq_csv("bar300K1ms_L5nm_xn")
if cx1:
    ax2.plot(cx1["u"], cx1["Tn"], "k-.", lw=1.6, label=r"$L=5$ nm (XN, 1 m/s)")
ax2.set_xlim(0, 0.72)
ax2.set_ylim(0, 17.0)
ax2.set_xlabel("Deslocamento total $u$ (nm)")
ax2.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax2.set_title(r"(b) Abaqus Riks a 300 K: $L_{\mathrm{snap}}^{\mathrm{XN}}=2.20$ nm vs. $L_{\mathrm{snap}}^{\mathrm{Bilin}}=3.89$ nm")
ax2.legend(loc="upper right", fontsize=7.5)
fig.tight_layout()
save_fig(fig, "fig02_abaqus_snapback_transition")

# FIG 3: LAMMPS vs Abaqus (0 K e 300 K)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.4))
cols_bar = {2: "C4", 5: "C0", 10: "C2", 20: "C1", 40: "C3"}
for L in [2, 5, 10, 20, 40]:
    md = np.loadtxt(BAR_DIR / f"bar0K_L{L}nm_full.dat", comments="#")
    ax1.plot(md[:, 0] * 0.1, md[:, 1], "-", color=cols_bar[L], lw=1.8, label=rf"LAMMPS $L={L}$ nm")
    abq = load_abq_csv(f"bar0K_L{L}nm_xn")
    if abq:
        ax1.plot(abq["u"], abq["Tn"], "--", color=cols_bar[L], lw=1.5, alpha=0.9, label=rf"Abaqus XN $L={L}$ nm")
ax1.set_xlim(0, 5.4)
ax1.set_ylim(-0.5, 17.5)
ax1.set_xlabel("Deslocamento total $u$ (nm)")
ax1.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax1.set_title("(a) Barras bicristal a 0 K: LAMMPS MS vs. Abaqus CZM (Riks)")
ax1.legend(ncol=2, fontsize=7.1, loc="lower right")

for L in [5, 10, 20, 40]:
    md = np.loadtxt(BAR_DIR / f"bar300K_L{L}nm_v10ms_full.dat", comments="#")
    tn_s = moving_avg(md[:, 1], 15)
    ax2.plot(md[:, 0] * 0.1, tn_s, "-", color=cols_bar[L], lw=1.7, label=rf"LAMMPS 10 m/s $L={L}$ nm")
    abq = load_abq_csv(f"bar300K_L{L}nm_xn")
    if abq:
        ax2.plot(abq["u"], abq["Tn"], "--", color=cols_bar[L], lw=1.5, alpha=0.9, label=rf"Abaqus XN $L={L}$ nm")
ax2.set_xlim(0, 6.2)
ax2.set_ylim(-2.5, 17.5)
ax2.set_xlabel("Deslocamento total $u$ (nm)")
ax2.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax2.set_title("(b) Barras bicristal a 300 K (10 m/s): LAMMPS MD vs. Abaqus CZM")
ax2.legend(ncol=2, fontsize=7.1, loc="lower right")
fig.tight_layout()
save_fig(fig, "fig03_lammps_vs_abaqus_bars")

# FIG 4: Efeito de escala u_pk(L) e taxa em L = 5 nm
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
L0_vals = np.array([2, 5, 10, 20, 40], dtype=float)
u0_md = np.array([bars0[f"L{int(L)}nm"]["u_peak_A"] * 0.1 for L in L0_vals])
u0_abq = np.array([load_abq_csv(f"bar0K_L{int(L)}nm_xn")["u"][np.argmax(load_abq_csv(f"bar0K_L{int(L)}nm_xn")["Tn"])] for L in L0_vals])

L3_vals = np.array([5, 10, 20, 40], dtype=float)
u3_md = np.array([bars300[f"L{int(L)}nm_v10ms"]["u_peak_A"] * 0.1 for L in L3_vals])
u3_abq = np.array([load_abq_csv(f"bar300K_L{int(L)}nm_xn")["u"][np.argmax(load_abq_csv(f"bar300K_L{int(L)}nm_xn")["Tn"])] for L in L3_vals])

ax1.plot(L0_vals, u0_md, "ko-", lw=1.8, ms=5.5, label=r"LAMMPS 0 K ($M_{\mathrm{sec}}=136.6$ GPa)")
ax1.plot(L0_vals, u0_abq, "k--s", lw=1.5, ms=5.0, mfc="w", label=r"Abaqus 0 K ($M_0=177.7$ GPa)")
ax1.plot(L3_vals, u3_md, "C3o-", lw=1.8, ms=5.5, label=r"LAMMPS 300 K, 10 m/s ($M_{\mathrm{sec}}=120.5$ GPa)")
ax1.plot(L3_vals, u3_abq, "C3--s", lw=1.5, ms=5.0, mfc="w", label=r"Abaqus 300 K ($M_0=171.2$ GPa)")
ax1.set_xlabel("Comprimento da barra $L$ (nm)")
ax1.set_ylabel(r"Deslocamento no pico $u_{\mathrm{pico}}$ (nm)")
ax1.set_title(r"(a) Escala de $u_{\mathrm{pico}}(L)$: efeito anarmónico nos grãos")
ax1.legend(loc="upper left", fontsize=7.6)

for tag, col, ls, lbl in [
    ("L5nm_v10ms",    "C0", "-",  "LAMMPS 300 K, 10 m/s (seed 1)"),
    ("L5nm_v10ms_s2", "C0", ":",  "LAMMPS 300 K, 10 m/s (seed 2)"),
    ("L5nm_v1ms",     "C1", "-",  "LAMMPS 300 K, 1 m/s"),
]:
    md = np.loadtxt(BAR_DIR / f"bar300K_{tag}_full.dat", comments="#")
    ax2.plot(md[:, 0]*0.1, moving_avg(md[:, 1], 15), ls, color=col, lw=1.7, label=lbl)

for jn, col, ls, lbl in [
    ("bar300K_L5nm_xn",    "C0", "--", "Abaqus XN 300 K (10 m/s)"),
    ("bar300K1ms_L5nm_xn", "C1", "--", "Abaqus XN 300 K (1 m/s)"),
]:
    abq = load_abq_csv(jn)
    if abq:
        ax2.plot(abq["u"], abq["Tn"], ls, color=col, lw=1.6, label=lbl)

ax2.set_xlim(0, 1.25)
ax2.set_ylim(-3.0, 17.2)
ax2.set_xlabel("Deslocamento total $u$ (nm)")
ax2.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax2.set_title(r"(b) Taxa de deformação e variabilidade estocástica ($L=5$ nm, 300 K)")
ax2.legend(loc="lower right", fontsize=7.5)
fig.tight_layout()
save_fig(fig, "fig04_scaling_and_rate")

# FIG 5: Convergencia de malha GCI (PBC vs Faces Livres)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.4, 4.3))
pbc_jobs = [("gci_pbc_M1", "M1 (66 el, h=0.50 nm)", "C0", ":"),
            ("bar0K_L5nm_xn", "M2 (420 el, h=0.25 nm)", "C0", "--"),
            ("gci_pbc_M3", "M3 (1488 el, h=0.17 nm)", "C0", "-")]
free_jobs = [("gci_free_M1", "M1 livre (66 el, h=0.50 nm)", "C3", ":"),
             ("gci_free_M2", "M2 livre (504 el, h=0.25 nm)", "C3", "--"),
             ("gci_free_M3", "M3 livre (1674 el, h=0.17 nm)", "C3", "-")]

for jn, lbl, col, ls in pbc_jobs + free_jobs:
    c = load_abq_csv(jn)
    if c:
        ax1.plot(c["u"], c["Tn"], ls, color=col, lw=1.7, label=lbl)
ax1.set_xlim(0, 1.05)
ax1.set_ylim(0, 17.0)
ax1.set_xlabel("Deslocamento total $u$ (nm)")
ax1.set_ylabel(r"Tração normal $T_n$ (GPa)")
ax1.set_title(r"(a) Convergência de malha e condição lateral ($L=5$ nm, 0 K)")
ax1.legend(loc="lower right", fontsize=7.3)

hs = np.array([0.50, 0.25, 1.0/6.0])
upk_pbc = [load_abq_csv(jn)["u"][np.argmax(load_abq_csv(jn)["Tn"])] for jn, _, _, _ in pbc_jobs]
upk_free = [load_abq_csv(jn)["u"][np.argmax(load_abq_csv(jn)["Tn"])] for jn, _, _, _ in free_jobs]

ax2.plot(hs, upk_pbc, "C0o-", lw=1.8, ms=6, label=r"PBC laterais ($\mathrm{GCI}_{21}=0.076\%$)")
ax2.plot(hs, upk_free, "C3s-", lw=1.8, ms=6, label=r"Faces livres ($\mathrm{GCI}_{21}=0.009\%$)")
ax2.invert_xaxis()
ax2.set_xlabel("Tamanho característico da malha $h$ (nm)")
ax2.set_ylabel(r"Deslocamento no pico $u_{\mathrm{pico}}$ (nm)")
ax2.set_title(r"(b) Refinamento de malha ($h = 0.50 \to 0.25 \to 0.17$ nm)")
ax2.legend(loc="center right", fontsize=8.0)
fig.tight_layout()
save_fig(fig, "fig05_gci_mesh_convergence")

# Tabela resumo
print("\n=== 0 K (LAMMPS MS vs Abaqus Xu-Needleman) ===")
print(f"{'L(nm)':>5} | {'t_MD(GPa)':>9} {'t_FE(GPa)':>9} {'Err_t':>7} | {'u_MD(nm)':>8} {'u_FE(nm)':>8} {'Err_u':>7} | {'Recuo_FE':>8}")
for L in [2, 5, 10, 20, 40]:
    t_md = bars0[f"L{L}nm"]["t_peak_GPa"]
    u_md = bars0[f"L{L}nm"]["u_peak_A"] * 0.1
    abq = load_abq_csv(f"bar0K_L{L}nm_xn")
    ip = int(np.argmax(abq["Tn"]))
    t_fe, u_fe = float(abq["Tn"][ip]), float(abq["u"][ip])
    rec = u_fe - float(np.min(abq["u"][ip:]))
    print(f"{L:5d} | {t_md:9.3f} {t_fe:9.3f} {(t_fe/t_md-1)*100:+6.2f}% | {u_md:8.4f} {u_fe:8.4f} {(u_fe/u_md-1)*100:+6.2f}% | {rec:8.4f}")

print("\n=== 300 K, 10 m/s (LAMMPS MD vs Abaqus Xu-Needleman) ===")
print(f"{'L(nm)':>5} | {'t_MD(GPa)':>9} {'t_FE(GPa)':>9} {'Err_t':>7} | {'u_MD(nm)':>8} {'u_FE(nm)':>8} {'Err_u':>7} | {'Recuo_FE':>8}")
for L in [5, 10, 20, 40]:
    t_md = bars300[f"L{L}nm_v10ms"]["t_peak_smooth_GPa"]
    u_md = bars300[f"L{L}nm_v10ms"]["u_peak_A"] * 0.1
    abq = load_abq_csv(f"bar300K_L{L}nm_xn")
    ip = int(np.argmax(abq["Tn"]))
    t_fe, u_fe = float(abq["Tn"][ip]), float(abq["u"][ip])
    rec = u_fe - float(np.min(abq["u"][ip:]))
    print(f"{L:5d} | {t_md:9.3f} {t_fe:9.3f} {(t_fe/t_md-1)*100:+6.2f}% | {u_md:8.4f} {u_fe:8.4f} {(u_fe/u_md-1)*100:+6.2f}% | {rec:8.4f}")

t_md1 = bars300["L5nm_v1ms"]["t_peak_smooth_GPa"]
u_md1 = bars300["L5nm_v1ms"]["u_peak_A"] * 0.1
abq1 = load_abq_csv("bar300K1ms_L5nm_xn")
ip1 = int(np.argmax(abq1["Tn"]))
print(f"\n=== 300 K, 1 m/s (L = 5 nm) ===\nt_MD={t_md1:.3f}, t_FE={abq1['Tn'][ip1]:.3f} GPa ({(abq1['Tn'][ip1]/t_md1-1)*100:+.2f}%) | u_MD={u_md1:.4f}, u_FE={abq1['u'][ip1]:.4f} nm ({(abq1['u'][ip1]/u_md1-1)*100:+.2f}%)")

print("\nFiguras geradas em 07_report_figs/:")
for p in sorted(FIG_DIR.glob("fig0*.pdf")):
    print("  -", p.name)

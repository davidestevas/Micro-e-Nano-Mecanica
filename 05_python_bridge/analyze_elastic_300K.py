import json
import numpy as np
from pathlib import Path

base = Path(__file__).resolve().parents[1] / "01_lammps_bulk"
delta = 0.01  # +/- 1.0%

def load_case(idx):
    d = np.loadtxt(base / f"stress_300K_case{idx}.dat", comments="#")
    # Cols: TimeStep sxx syy szz syz sxz sxy T (dt_bin = 0.01 ps)
    return d[d[:, 0] > 0]

def block_se(y, dt_bin=0.01, block_counts=(100, 50, 25, 10, 5)):
    info = []
    for nb in block_counts:
        lb = len(y) // nb
        bmeans = y[:nb * lb].reshape(nb, lb).mean(axis=1)
        se = float(np.std(bmeans, ddof=1) / np.sqrt(nb))
        info.append((nb, lb * dt_bin, se))
    # Patamar conservador entre blocos de 1 a 10 ps
    se_plateau = max(x[2] for x in info if x[1] >= 1.0)
    return float(np.mean(y)), se_plateau, info

c1, c2, c3, c4 = [load_case(i) for i in (1, 2, 3, 4)]

print("=== Temperaturas medias nos 4 casos NVT (50 ps) ===")
for i, c in enumerate([c1, c2, c3, c4], start=1):
    print(f"  Caso {i}: <T> = {np.mean(c[:, 7]):.2f} +/- {np.std(c[:, 7], ddof=1):.2f} K")

def central_diff_const(y_pos, y_neg, name):
    mp, sep, infop = block_se(y_pos)
    mm, sem, infom = block_se(y_neg)
    val = (mp - mm) / (2.0 * delta)
    se = np.sqrt(sep**2 + sem**2) / (2.0 * delta)
    print(f"\n--- {name} ---")
    print(f"  <sigma+> = {mp:+.4f} +/- {sep:.4f} GPa | <sigma-> = {mm:+.4f} +/- {sem:.4f} GPa")
    print(f"  {name} = {val:.2f} GPa  (SE = {se:.2f} GPa, 95% CI: [{val - 1.96*se:.2f}, {val + 1.96*se:.2f}] GPa)")
    for (nb, lb_ps, sp), (_, _, sm) in zip(infop, infom):
        se_b = np.sqrt(sp**2 + sm**2) / (2.0 * delta)
        print(f"    N_b={nb:3d} (l_b={lb_ps:5.2f} ps): SE({name}) = {se_b:.3f} GPa")
    return val, se

C11, se11 = central_diff_const(c1[:, 1], c2[:, 1], "C11")
C12, se12 = central_diff_const(c1[:, 2], c2[:, 2], "C12")
C13, se13 = central_diff_const(c1[:, 3], c2[:, 3], "C13")
C44, se44 = central_diff_const(c3[:, 6], c4[:, 6], "C44")

C12_avg = 0.5 * (C12 + C13)
se12_avg = 0.5 * np.sqrt(se12**2 + se13**2)

# Verificacao de estabilidade de Born, modulo volumetrico K e anisotropia de Zener A
born_ok = (C11 - C12_avg > 0) and (C44 > 0) and (C11 + 2 * C12_avg > 0)
K = (C11 + 2.0 * C12_avg) / 3.0
A = 2.0 * C44 / (C11 - C12_avg)

# Rotacao do tensor cubico para a orientacao [310] (normal z = [310]/sqrt(10))
def rotate_310(c11, c12, c44):
    C = np.zeros((3, 3, 3, 3))
    for i in range(3):
        C[i, i, i, i] = c11
        for j in range(3):
            if i != j:
                C[i, i, j, j] = c12
                C[i, j, i, j] = c44
                C[i, j, j, i] = c44
    e1 = np.array([-1.0, 3.0, 0.0]) / np.sqrt(10.0)
    e2 = np.array([0.0, 0.0, 1.0])
    e3 = np.array([3.0, 1.0, 0.0]) / np.sqrt(10.0)
    R = np.vstack([e1, e2, e3])
    Cr = np.einsum("ia,jb,kc,ld,abcd->ijkl", R, R, R, R, C)
    vidx = [(0,0), (1,1), (2,2), (1,2), (0,2), (0,1)]
    CV = np.array([[Cr[i, j, k, l] for (k, l) in vidx] for (i, j) in vidx])
    SV = np.linalg.inv(CV)
    C33 = CV[2, 2]
    C35 = CV[2, 4]
    C55 = CV[4, 4]
    M_free_shear = C33 - C35**2 / C55
    E_310 = 1.0 / SV[2, 2]
    return C33, M_free_shear, E_310, CV

C33_0K, Mfs_0K, E310_0K, _ = rotate_310(169.92, 122.63, 76.21)
C33_300K, Mfs_300K, E310_300K, _ = rotate_310(C11, C12_avg, C44)

print("\n============================================================")
print("RESUMO ELASTICO DE COBRE (Mishin 2001): 0 K vs 300 K")
print("============================================================")
print(f"Simetria cubica |C12 - C13| = {abs(C12 - C13):.2f} GPa")
print(f"Condicoes de Born           = {'OK' if born_ok else 'FAIL'}")
print(f"0 K   : C11=169.92, C12=122.63, C44=76.21 GPa | K=138.39 GPa, A=3.223")
print(f"300 K : C11={C11:.2f} +/- {se11:.2f}, C12={C12_avg:.2f} +/- {se12_avg:.2f}, C44={C44:.2f} +/- {se44:.2f} GPa | K={K:.2f} GPa, A={A:.3f}")
print(f"Modulos rodados na normal [310]:")
print(f"  0 K   -> C'_3333 = {C33_0K:.2f} GPa | M_free_shear = {Mfs_0K:.2f} GPa | E_[310] = {E310_0K:.2f} GPa")
print(f"  300 K -> C'_3333 = {C33_300K:.2f} GPa | M_free_shear = {Mfs_300K:.2f} GPa | E_[310] = {E310_300K:.2f} GPa")

out = {
    "T_0K": {
        "a0_A": 3.614925, "Ecoh_eV": -3.540218, "Evf_eV": 1.2726,
        "C11_GPa": 169.92, "C12_GPa": 122.63, "C44_GPa": 76.21,
        "C3333_310_GPa": C33_0K, "M_free_shear_310_GPa": Mfs_0K, "E_310_GPa": E310_0K
    },
    "T_300K": {
        "a_A": 3.631868, "a_SE_A": 0.000024,
        "C11_GPa": C11, "C11_SE_GPa": se11,
        "C12_GPa": C12_avg, "C12_SE_GPa": se12_avg,
        "C13_check_GPa": C13,
        "C44_GPa": C44, "C44_SE_GPa": se44,
        "K_GPa": K, "A_zener": A,
        "C3333_310_GPa": C33_300K, "M_free_shear_310_GPa": Mfs_300K, "E_310_GPa": E310_300K
    }
}
out_path = Path(__file__).resolve().parent / "cu_elastic_summary.json"
out_path.write_text(json.dumps(out, indent=2))
print(f"\nGuardado em: {out_path}")

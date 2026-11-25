from __future__ import print_function, division
import math

def make_c1_bilinear_table(t0, Gc, Kn, frac_fillet=0.10):
    d0_nom = t0 / Kn
    w = frac_fillet * d0_nom
    Ks = t0 / (2.0 * Gc / t0 - d0_nom)
    for _ in range(50):
        d1 = (t0 - 0.5 * (Kn ** 2) / (Kn + Ks) * w) / Kn
        d2 = d1 + w
        T2 = Kn * d1 + 0.5 * (Kn - Ks) * w
        df = d2 + T2 / Ks
        A1 = 0.5 * Kn * d1 * d1
        A2 = Kn * d1 * w + 0.5 * Kn * w * w - (Kn + Ks) * w * w / 6.0
        A3 = 0.5 * T2 * T2 / Ks
        Ks *= ((A1 + A2 + A3) / Gc)

    t_init = Kn * d1
    rows = ["0.00000000, 0.00000000"]
    for i in range(1, 41):
        d = d1 + (i / 40.0) * w
        T = Kn * d1 + Kn * (d - d1) - 0.5 * (Kn + Ks) / w * ((d - d1) ** 2)
        D = 1.0 - T / (Kn * d)
        rows.append("%.8f, %.8f" % (D, d - d1))
    for i in range(1, 81):
        d = d2 + (i / 80.0) * (df - d2)
        T = max(0.0, T2 - Ks * (d - d2))
        D = (1.0 - T / (Kn * d)) if i < 80 else 0.99990000
        rows.append("%.8f, %.8f" % (D, d - d1))
    return t_init, "\n".join(rows)

def build_z_half(L_half, h_fine=0.25, z_fine=2.0, h_coarse=1.0):
    if L_half <= z_fine + 1.0e-6:
        n = max(2, int(round(L_half / h_fine)))
        return [L_half * i / n for i in range(n + 1)]
    n1 = max(2, int(round(z_fine / h_fine)))
    z1 = [z_fine * i / n1 for i in range(n1 + 1)]
    rem = L_half - z_fine
    n2 = max(1, int(round(rem / h_coarse)))
    z2 = [z_fine + rem * i / n2 for i in range(1, n2 + 1)]
    return z1 + z2

def write_bar_inp(jobname, L_nm, law="bilin", temp_K=0, nx=4, ny=5, h_fine=0.25, z_fine=2.0, h_coarse=1.0, bc="pbc", riks_inc=300, arc_max=0.020, frac_fillet=0.10):
    if temp_K == 0:
        Lx, Ly = 1.143139657, 1.445970000
        C11, C12, C44 = 169.92, 122.63, 76.21
        t0, Gc, Kn_bil = 15.6991, 2.126821, 545.70
    elif temp_K == 30010:
        Lx, Ly = 1.148498, 1.452747
        C11, C12, C44 = 163.89, 120.15, 71.68
        t0, Gc, Kn_bil = 15.909, 3.2449, 340.40
    elif temp_K == 30001:
        Lx, Ly = 1.148498, 1.452747
        C11, C12, C44 = 163.89, 120.15, 71.68
        t0, Gc, Kn_bil = 15.025, 3.0626, 390.80
    else:
        raise ValueError("temp_K invalido")

    Area = Lx * Ly
    F_ref = t0 * Area

    z_pos = build_z_half(0.5 * L_nm, h_fine=h_fine, z_fine=z_fine, h_coarse=h_coarse)
    nz_half = len(z_pos) - 1
    z_g1 = [-z for z in reversed(z_pos)]
    z_g2 = list(z_pos)
    z_all = z_g1 + z_g2
    nk = len(z_all)

    dx = Lx / nx
    dy = Ly / ny

    def nid(i, j, k):
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1

    rp_id = nk * (nx + 1) * (ny + 1) + 1

    out = []
    out.append("*HEADING")
    out.append("Cu Sigma5(310) 3D bar | L=%.3f nm | law=%s | T=%d | bc=%s | mesh=%dx%dx%d | units=nm-nN-GPa" % (
        L_nm, law, temp_K, bc, nx, ny, 2 * nz_half))

    out.append("*NODE")
    for k in range(nk):
        z = z_all[k]
        for j in range(ny + 1):
            y = j * dy
            for i in range(nx + 1):
                out.append("%d, %.10f, %.10f, %.10f" % (nid(i, j, k), i * dx, y, z))
    out.append("%d, %.10f, %.10f, %.10f" % (rp_id, 0.5 * Lx, 0.5 * Ly, 0.5 * L_nm + 0.5))

    eid = 0
    out.append("*ELEMENT, TYPE=C3D8, ELSET=GRAIN1")
    for k in range(nz_half):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                out.append("%d, %d, %d, %d, %d, %d, %d, %d, %d" % (
                    eid, nid(i,j,k), nid(i+1,j,k), nid(i+1,j+1,k), nid(i,j+1,k),
                    nid(i,j,k+1), nid(i+1,j,k+1), nid(i+1,j+1,k+1), nid(i,j+1,k+1)))

    out.append("*ELEMENT, TYPE=C3D8, ELSET=GRAIN2")
    for k in range(nz_half + 1, 2 * nz_half + 1):
        for j in range(ny):
            for i in range(nx):
                eid += 1
                out.append("%d, %d, %d, %d, %d, %d, %d, %d, %d" % (
                    eid, nid(i,j,k), nid(i+1,j,k), nid(i+1,j+1,k), nid(i,j+1,k),
                    nid(i,j,k+1), nid(i+1,j,k+1), nid(i+1,j+1,k+1), nid(i,j+1,k+1)))

    out.append("*ELEMENT, TYPE=COH3D8, ELSET=COH")
    k_lo, k_hi = nz_half, nz_half + 1
    for j in range(ny):
        for i in range(nx):
            eid += 1
            out.append("%d, %d, %d, %d, %d, %d, %d, %d, %d" % (
                eid, nid(i,j,k_lo), nid(i+1,j,k_lo), nid(i+1,j+1,k_lo), nid(i,j+1,k_lo),
                nid(i,j,k_hi), nid(i+1,j,k_hi), nid(i+1,j+1,k_hi), nid(i,j+1,k_hi)))

    def append_nset(name, node_list):
        out.append("*NSET, NSET=%s" % name)
        for p in range(0, len(node_list), 16):
            out.append(", ".join(str(x) for x in node_list[p:p + 16]))

    bot_nodes = [nid(i, j, 0) for j in range(ny + 1) for i in range(nx + 1)]
    top_nodes = [nid(i, j, nk - 1) for j in range(ny + 1) for i in range(nx + 1)]
    ic, jc = nx // 2, ny // 2
    append_nset("BOT", bot_nodes)
    append_nset("TOP", top_nodes)
    append_nset("RP", [rp_id])
    append_nset("IF_BOT_0", [nid(ic, jc, k_lo)])
    append_nset("IF_TOP_0", [nid(ic, jc, k_hi)])

    if bc == "pbc":
        for k in range(1, nk - 1):
            for dof in (1, 2, 3):
                out.append("*EQUATION\n2\n%d, %d, 1.0, %d, %d, -1.0" % (nid(nx, ny, k), dof, nid(0, 0, k), dof))
            for j in range(ny):
                for dof in (1, 2, 3):
                    out.append("*EQUATION\n2\n%d, %d, 1.0, %d, %d, -1.0" % (nid(nx, j, k), dof, nid(0, j, k), dof))
            for i in range(nx):
                for dof in (1, 2, 3):
                    out.append("*EQUATION\n2\n%d, %d, 1.0, %d, %d, -1.0" % (nid(i, ny, k), dof, nid(i, 0, k), dof))

    for nd in top_nodes:
        out.append("*EQUATION\n2\n%d, 3, 1.0, %d, 3, -1.0" % (nd, rp_id))

    s10 = math.sqrt(10.0)
    out.append("*ORIENTATION, NAME=ORI_G1, DEFINITION=COORDINATES")
    out.append("%.12f, 0.0, %.12f, %.12f, 0.0, %.12f" % (-1.0/s10, 3.0/s10, 3.0/s10, 1.0/s10))
    out.append("3, 0.0")
    out.append("*ORIENTATION, NAME=ORI_G2, DEFINITION=COORDINATES")
    out.append("%.12f, 0.0, %.12f, %.12f, 0.0, %.12f" % ( 1.0/s10, 3.0/s10, 3.0/s10, -1.0/s10))
    out.append("3, 0.0")

    out.append("*SOLID SECTION, ELSET=GRAIN1, MATERIAL=CU_BULK, ORIENTATION=ORI_G1")
    out.append("*SOLID SECTION, ELSET=GRAIN2, MATERIAL=CU_BULK, ORIENTATION=ORI_G2")
    out.append("*COHESIVE SECTION, ELSET=COH, MATERIAL=MAT_COH, RESPONSE=TRACTION SEPARATION, THICKNESS=SPECIFIED")
    out.append("1.0")

    out.append("*MATERIAL, NAME=CU_BULK")
    out.append("*ELASTIC, TYPE=ORTHOTROPIC")
    out.append("%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f" % (C11, C12, C11, C12, C12, C11, C44, C44))
    out.append("%.4f" % C44)

    out.append("*MATERIAL, NAME=MAT_COH")
    if law == "bilin":
        t_init, tab_str = make_c1_bilinear_table(t0, Gc, Kn_bil, frac_fillet=frac_fillet)
        out.append("*ELASTIC, TYPE=TRACTION")
        out.append("%.4f, %.4f, %.4f" % (Kn_bil, Kn_bil, Kn_bil))
        out.append("*DAMAGE INITIATION, CRITERION=MAXS")
        out.append("%.6f, %.6f, %.6f" % (t_init, t_init, t_init))
        out.append("*DAMAGE EVOLUTION, TYPE=DISPLACEMENT, SOFTENING=TABULAR")
        out.append(tab_str)
    elif law == "xn":
        dc = Gc / (math.e * t0)
        Kn_xn = math.e * t0 / dc
        d_init = 0.005 * dc
        t_init = Kn_xn * d_init
        out.append("*ELASTIC, TYPE=TRACTION")
        out.append("%.4f, %.4f, %.4f" % (Kn_xn, Kn_xn, Kn_xn))
        out.append("*DAMAGE INITIATION, CRITERION=MAXS")
        out.append("%.6f, %.6f, %.6f" % (t_init, t_init, t_init))
        out.append("*DAMAGE EVOLUTION, TYPE=DISPLACEMENT, SOFTENING=TABULAR")
        out.append("0.00000000, 0.00000000")
        for i in range(1, 121):
            s = (i / 120.0) ** 1.6
            delta = d_init + s * (10.0 * dc - d_init)
            D_val = (1.0 - math.exp(-delta / dc)) if i < 120 else 0.99990000
            out.append("%.8f, %.8f" % (D_val, delta - d_init))
    else:
        raise ValueError("law invalido")

    out.append("*STEP, NAME=PULL, NLGEOM=NO, INC=%d" % riks_inc)
    out.append("*STATIC, RIKS")
    out.append("0.01, 1.0, 1e-8, %.4f" % arc_max)
    out.append("*BOUNDARY")
    out.append("BOT, 1, 3, 0.0")
    out.append("TOP, 1, 2, 0.0")
    out.append("RP, 1, 2, 0.0")
    out.append("*CLOAD")
    out.append("RP, 3, %.6f" % F_ref)
    out.append("*OUTPUT, FIELD, FREQUENCY=25")
    out.append("*NODE OUTPUT")
    out.append("U, RF")
    out.append("*ELEMENT OUTPUT, ELSET=COH")
    out.append("S, SDEG")
    out.append("*OUTPUT, HISTORY, FREQUENCY=1")
    out.append("*NODE OUTPUT, NSET=RP")
    out.append("U3,")
    out.append("*NODE OUTPUT, NSET=IF_BOT_0")
    out.append("U3,")
    out.append("*NODE OUTPUT, NSET=IF_TOP_0")
    out.append("U3,")
    out.append("*ENERGY OUTPUT")
    out.append("ALLWK, ALLSE, ALLDMD, ALLVD")
    out.append("*END STEP")

    fname = jobname + ".inp"
    with open(fname, "w") as f:
        f.write("\n".join(out) + "\n")
    print("Escrito %-26s | L=%5.1f nm | el=%4d | bc=%s" % (fname, L_nm, eid, bc))

if __name__ == "__main__":
    # 1) 3 correcoes da Bateria 1 (0 K)
    write_bar_inp("bar0K_L10nm_bilin", 10.0, law="bilin", temp_K=0, riks_inc=350, arc_max=0.015, frac_fillet=0.15)
    write_bar_inp("bar0K_L40nm_bilin", 40.0, law="bilin", temp_K=0, riks_inc=350, arc_max=0.015, frac_fillet=0.15)
    write_bar_inp("bar0K_L20nm_xn",    20.0, law="xn",    temp_K=0, riks_inc=350, arc_max=0.010)

    # 2) Barras a 300 K (10 m/s: L=2,3,5,10,20,40 nm XN + L=2,3,5 nm Bilin; 1 m/s: L=5 nm XN & Bilin)
    for L_nm in [2.0, 3.0, 5.0, 10.0, 20.0, 40.0]:
        amax = 0.012 if L_nm >= 20.0 else 0.020
        write_bar_inp("bar300K_L%dnm_xn" % int(L_nm), L_nm, law="xn", temp_K=30010, riks_inc=320, arc_max=amax)
    for L_nm in [2.0, 3.0, 5.0]:
        write_bar_inp("bar300K_L%dnm_bilin" % int(L_nm), L_nm, law="bilin", temp_K=30010, riks_inc=300, arc_max=0.015, frac_fillet=0.15)
    write_bar_inp("bar300K1ms_L5nm_xn",    5.0, law="xn",    temp_K=30001, riks_inc=300, arc_max=0.020)
    write_bar_inp("bar300K1ms_L5nm_bilin", 5.0, law="bilin", temp_K=30001, riks_inc=300, arc_max=0.015, frac_fillet=0.15)

    # 3) Estudo de convergencia em 3 malhas (GCI) em L=5 nm, 0 K, Xu-Needleman (faces livres + check PBC)
    write_bar_inp("gci_free_M1", 5.0, law="xn", temp_K=0, nx=2, ny=3, h_fine=0.500, z_fine=2.5, h_coarse=0.500, bc="free", riks_inc=300, arc_max=0.015)
    write_bar_inp("gci_free_M2", 5.0, law="xn", temp_K=0, nx=4, ny=6, h_fine=0.250, z_fine=2.5, h_coarse=0.250, bc="free", riks_inc=300, arc_max=0.015)
    write_bar_inp("gci_free_M3", 5.0, law="xn", temp_K=0, nx=6, ny=9, h_fine=0.167, z_fine=2.5, h_coarse=0.167, bc="free", riks_inc=300, arc_max=0.015)
    write_bar_inp("gci_pbc_M1",  5.0, law="xn", temp_K=0, nx=2, ny=3, h_fine=0.500, z_fine=2.5, h_coarse=0.500, bc="pbc",  riks_inc=300, arc_max=0.015)
    write_bar_inp("gci_pbc_M3",  5.0, law="xn", temp_K=0, nx=6, ny=8, h_fine=0.167, z_fine=2.5, h_coarse=0.167, bc="pbc",  riks_inc=300, arc_max=0.015)

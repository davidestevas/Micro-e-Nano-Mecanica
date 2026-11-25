import math

a0 = 3.614925
px = math.sqrt(10.0) / 2.0 * a0   # 5.7156982853 A
py = a0                           # 3.6149250000 A
pz = math.sqrt(10.0) / 2.0 * a0   # 5.7156982853 A
dz = pz / 10.0                    # 0.5715698285 A

nx, ny = 2, 4
Lx, Ly = nx * px, ny * py

def build_bar(L_nm):
    L_A = L_nm * 10.0
    # Metade do comprimento movel + 12.5 A de garra rigida (> 2*r_cut = 11.01 A)
    nz_grain = int(math.ceil((L_A / 2.0 + 12.5) / pz))
    n_planes = 10 * nz_grain
    atoms = []
    for k in range(n_planes + 1):
        ux = ((-3 * k) % 10) / 10.0
        uy = 0.5 * (k % 2)
        for ix in range(nx):
            for iy in range(ny):
                x = (ix + ux) * px
                y = (iy + uy) * py
                atoms.append((1, x, y, -k * dz))
                atoms.append((2, x, y, +k * dz))
    H = n_planes * dz
    zlo, zhi = -H - 15.0, +H + 15.0
    fname = f"data.bar_L{L_nm}nm_init"
    with open(fname, "w") as f:
        f.write(f"Cu Sigma5(310)[001] bicrystal bar L={L_nm} nm (H_grain={H:.2f} A)\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write("2 atom types\n\n")
        f.write(f"0.0000000000 {Lx:.10f} xlo xhi\n")
        f.write(f"0.0000000000 {Ly:.10f} ylo yhi\n")
        f.write(f"{zlo:.10f} {zhi:.10f} zlo zhi\n\n")
        f.write("Masses\n\n")
        f.write("1 63.546000 # Cu\n")
        f.write("2 63.546000 # Cu\n\n")
        f.write("Atoms # atomic\n\n")
        for idx, (atype, x, y, z) in enumerate(atoms, start=1):
            f.write(f"{idx} {atype} {x:.10f} {y:.10f} {z:.10f}\n")
    print(f"Escrito {fname}: {len(atoms):4d} atomos | H_grain={H:6.2f} A | zmob={L_A/2.0:5.1f} A")

if __name__ == "__main__":
    for L_nm in [2, 5, 10, 20, 40]:
        build_bar(L_nm)

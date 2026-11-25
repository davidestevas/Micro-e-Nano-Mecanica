import math
from pathlib import Path

a0 = 3.614925
px = math.sqrt(10.0) / 2.0 * a0   # 5.7156982853 A (periodo primitivo [-1 3 0])
py = a0                           # 3.6149250000 A (eixo de tilt [0 0 1])
pz = math.sqrt(10.0) / 2.0 * a0   # 5.7156982853 A (periodo normal [3 1 0])
dz = pz / 10.0                    # 0.5715698285 A (espacamento entre planos (310))

def write_lammps_data(fname, atoms, Lx, Ly, zlo, zhi, ntypes, title):
    with open(fname, "w") as f:
        f.write(f"{title}\n\n")
        f.write(f"{len(atoms)} atoms\n")
        f.write(f"{ntypes} atom types\n\n")
        f.write(f"0.0000000000 {Lx:.10f} xlo xhi\n")
        f.write(f"0.0000000000 {Ly:.10f} ylo yhi\n")
        f.write(f"{zlo:.10f} {zhi:.10f} zlo zhi\n\n")
        f.write("Masses\n\n")
        for t in range(1, ntypes + 1):
            f.write(f"{t} 63.546000 # Cu\n")
        f.write("\nAtoms # atomic\n\n")
        for idx, (atype, x, y, z) in enumerate(atoms, start=1):
            f.write(f"{idx} {atype} {x:.10f} {y:.10f} {z:.10f}\n")
    print(f"Escrito {fname}: {len(atoms)} atomos | Lx={Lx:.4f}, Ly={Ly:.4f}, z=[{zlo:.2f}, {zhi:.2f}] A")

def build_surface_slab(nx, ny, nz_cells, fname):
    # Placa com duas superficies livres (310) e 20 A de vacuo de cada lado
    Lx, Ly = nx * px, ny * py
    n_planes = 10 * nz_cells
    atoms = []
    for k in range(n_planes):
        ux = ((3 * k) % 10) / 10.0
        uy = 0.5 * (k % 2)
        z = k * dz
        for ix in range(nx):
            for iy in range(ny):
                atoms.append((1, (ix + ux) * px, (iy + uy) * py, z))
    zlo = -20.0
    zhi = (n_planes - 1) * dz + 20.0
    write_lammps_data(fname, atoms, Lx, Ly, zlo, zhi, 1, f"Cu (310) free surface slab ({n_planes} planes)")

def build_bicrystal(nx, ny, nz_grain, fname):
    # Bicristal simetrico Sigma5(310)[001]:
    # Tipo 1 (grao inferior): z <= 0 (k = 0..n_planes em z = -k*dz)
    # Tipo 2 (grao superior, espelho exato em z=0): z >= 0 (k = 0..n_planes em z = +k*dz)
    Lx, Ly = nx * px, ny * py
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
    zlo = -H - 15.0
    zhi = +H + 15.0
    write_lammps_data(fname, atoms, Lx, Ly, zlo, zhi, 2, f"Cu Sigma5(310)[001] bicrystal H_grain={H:.2f} A")

if __name__ == "__main__":
    # 1. Superficies (310): 40 planos (~22.9 A) e 80 planos (~45.7 A)
    build_surface_slab(2, 4, 4, "data.surf310_40p")
    build_surface_slab(2, 4, 8, "data.surf310_80p")

    # 2. Bicristais Sigma5(310): H_grain = 40.0 A (7*pz) e H_grain = 80.0 A (14*pz, dobro)
    build_bicrystal(2, 4, 7, "data.gb_H40")
    build_bicrystal(2, 4, 14, "data.gb_H80")

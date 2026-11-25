import numpy as np
from pathlib import Path

for Lnm, u_start in [(20, 24.0), (40, 42.0)]:
    dump_path = Path(f"dump.bar0K_L{Lnm}nm.lammpstrj")
    lines = dump_path.read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    blk = lines[starts[-1]:]

    n_at = int(blk[3].strip())
    xlo, xhi = map(float, blk[5].split()[:2])
    ylo, yhi = map(float, blk[6].split()[:2])
    zlo, zhi = map(float, blk[7].split()[:2])

    atoms = []
    for r in blk[9:9 + n_at]:
        p = r.split()
        atoms.append((int(p[0]), int(p[1]), float(p[2]), float(p[3]), float(p[4])))
    atoms.sort(key=lambda t: t[0])

    out_name = f"data.bar_L{Lnm}nm_cont_init"
    with open(out_name, "w") as f:
        f.write(f"Continuation state for L={Lnm} nm at u_start={u_start:.2f} A\n\n")
        f.write(f"{n_at} atoms\n2 atom types\n\n")
        f.write(f"{xlo:.10f} {xhi:.10f} xlo xhi\n")
        f.write(f"{ylo:.10f} {yhi:.10f} ylo yhi\n")
        f.write(f"{zlo - 15.0:.10f} {zhi + 30.0:.10f} zlo zhi\n\n")
        f.write("Masses\n\n1 63.546000 # Cu\n2 63.546000 # Cu\n\n")
        f.write("Atoms # atomic\n\n")
        for aid, atype, x, y, z in atoms:
            f.write(f"{aid} {atype} {x:.10f} {y:.10f} {z:.10f}\n")

    # Ler ultimo dE_Jm2 do ficheiro .dat original para manter continuidade exata
    arr = np.loadtxt(f"bar0K_L{Lnm}nm.dat", comments="#")
    u_last, tn_last, dE_last = arr[-1, 1], arr[-1, 2], arr[-1, 5]
    print(f"L={Lnm}nm -> {out_name} ({n_at} atomos) | u_last={u_last:.2f} A, tn_last={tn_last:.4f} GPa, dE_last={dE_last:.6f} J/m^2")

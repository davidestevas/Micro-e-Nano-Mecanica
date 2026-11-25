import numpy as np
from pathlib import Path

d_dir = Path(__file__).resolve().parents[1] / "03_lammps_czm_pull"

for vtag in ["v10ms", "v1ms"]:
    print(f"\n================ {vtag} ================")
    print((d_dir / f"pull300K_{vtag}_meta.txt").read_text().strip())
    arr = np.loadtxt(d_dir / f"pull300K_{vtag}.dat", comments="#")
    u, tn, tsx, tsy = arr[:, 1], arr[:, 2], arr[:, 5], arr[:, 6]
    print("Amostragem T_n(u) ao longo do ensaio:")
    for u_target in [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
        idx = int(np.argmin(np.abs(u - u_target)))
        w = slice(max(0, idx - 5), min(len(u), idx + 6))
        print(f"  u={u[idx]:4.2f} A -> Tn={np.mean(tn[w]):+7.3f} GPa | tsx={np.mean(tsx[w]):+6.3f} GPa")

    # Ler primeiro e ultimo frame do dump
    lines = (d_dir / f"dump.pull300K_{vtag}.lammpstrj").read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]
    for lbl, st in [("Primeiro frame", starts[0]), ("Ultimo frame", starts[-1])]:
        end = starts[starts.index(st) + 1] if starts.index(st) + 1 < len(starts) else len(lines)
        blk = lines[st:end]
        ai = [i for i, l in enumerate(blk) if l.startswith("ITEM: ATOMS")][0]
        pts = np.array([[float(x) for x in r.split()[:5]] for r in blk[ai + 1:] if r.strip()])
        t1 = pts[pts[:, 1] == 1, 4]
        t2 = pts[pts[:, 1] == 2, 4]
        zs = np.sort(pts[:, 4])
        gaps = np.diff(zs)
        ig = int(np.argmax(gaps))
        print(f"  {lbl}: z_all=[{zs[0]:.2f}, {zs[-1]:.2f}] A | tipo1_max={t1.max():.2f}, tipo2_min={t2.min():.2f} A | maior gap={gaps[ig]:.2f} A em z=[{zs[ig]:.2f}, {zs[ig+1]:.2f}] A")

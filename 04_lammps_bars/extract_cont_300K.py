import re
import numpy as np
from pathlib import Path

def parse_frame(lines, start_idx):
    n_at = int(lines[start_idx + 3].strip())
    xlo, xhi = map(float, lines[start_idx + 5].split()[:2])
    ylo, yhi = map(float, lines[start_idx + 6].split()[:2])
    zlo, zhi = map(float, lines[start_idx + 7].split()[:2])
    step = int(lines[start_idx + 1].strip())
    data = np.array([[float(x) for x in r.split()[:5]] for r in lines[start_idx + 9:start_idx + 9 + n_at]])
    data = data[np.argsort(data[:, 0])]
    return step, data, (xlo, xhi, ylo, yhi, zlo, zhi)

for Lnm, frame_offset in [(20, -3), (40, -1)]:
    tag = f"L{Lnm}nm_v10ms"
    lines = Path(f"dump.bar300K_{tag}.lammpstrj").read_text().splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("ITEM: TIMESTEP")]

    step0, d0, _ = parse_frame(lines, starts[0])
    step_c, dc, box = parse_frame(lines, starts[frame_offset])
    u0 = step_c * 0.001 * 0.10  # dt=0.001 ps, vel=0.10 A/ps

    dz = dc[:, 4] - d0[:, 4]
    bot_mask = np.abs(dz) < 1e-3
    top_mask = np.abs(dz - u0) < 1e-3
    z_bot_cut = float(dc[bot_mask, 4].max()) + 0.05
    z_top_cut = float(dc[top_mask, 4].min()) - 0.05

    data_out = f"data.bar300K_L{Lnm}nm_cont_init"
    xlo, xhi, ylo, yhi, zlo, zhi = box
    with open(data_out, "w") as f:
        f.write(f"300 K continuation L={Lnm} nm at u0={u0:.4f} A\n\n")
        f.write(f"{len(dc)} atoms\n2 atom types\n\n")
        f.write(f"{xlo:.10f} {xhi:.10f} xlo xhi\n")
        f.write(f"{ylo:.10f} {yhi:.10f} ylo yhi\n")
        f.write(f"{zlo - 15.0:.10f} {zhi + 30.0:.10f} zlo zhi\n\n")
        f.write("Masses\n\n1 63.546000 # Cu\n2 63.546000 # Cu\n\n")
        f.write("Atoms # atomic\n\n")
        for row in dc:
            f.write(f"{int(row[0])} {int(row[1])} {row[2]:.10f} {row[3]:.10f} {row[4]:.10f}\n")

    inc_out = f"cont_300K_L{Lnm}nm.inc"
    Path(inc_out).write_text(
        f"variable u0 equal {u0:.6f}\n"
        f"variable z_bot equal {z_bot_cut:.6f}\n"
        f"variable z_top equal {z_top_cut:.6f}\n"
    )
    print(f"L={Lnm}nm: step={step_c} (u0={u0:.2f} A) | Nbot={bot_mask.sum()}, Ntop={top_mask.sum()}, Nmob={(~bot_mask & ~top_mask).sum()}")

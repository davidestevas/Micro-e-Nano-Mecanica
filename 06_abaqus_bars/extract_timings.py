from __future__ import print_function
import glob, os, re

rows = []
for dat in sorted(glob.glob("*.dat")):
    job = os.path.splitext(dat)[0]
    txt = open(dat, "r").read()
    m_cpu = re.search(r"TOTAL CPU TIME \(SEC\)\s*=\s*([0-9.]+)", txt)
    m_wall = re.search(r"WALLCLOCK TIME \(SEC\)\s*=\s*([0-9.]+)", txt)
    cpu = float(m_cpu.group(1)) if m_cpu else -1.0
    wall = float(m_wall.group(1)) if m_wall else -1.0
    sta_ok = "UNKNOWN"
    if os.path.exists(job + ".sta"):
        st = open(job + ".sta", "r").read()
        sta_ok = "COMPLETED" if "COMPLETED SUCCESSFULLY" in st else "INCOMPLETE"
    rows.append((job, sta_ok, cpu, wall))

with open("abaqus_timings.csv", "w") as f:
    f.write("job,status,cpu_sec,wallclock_sec\n")
    for r in rows:
        f.write("%s,%s,%.2f,%.2f\n" % r)
print("Gravado abaqus_timings.csv com %d jobs." % len(rows))

from __future__ import print_function, division
import sys
import os
from odbAccess import openOdb

def job_params(jobname):
    if "300K1ms" in jobname:
        Lx, Ly, t0 = 1.148498, 1.452747, 15.025
    elif "300K" in jobname:
        Lx, Ly, t0 = 1.148498, 1.452747, 15.909
    else:
        Lx, Ly, t0 = 1.143139657, 1.445970000, 15.6991
    Area = Lx * Ly
    return Area, t0 * Area

def first_node_label(odb, nset_name):
    if nset_name in odb.rootAssembly.nodeSets:
        ns = odb.rootAssembly.nodeSets[nset_name].nodes[0]
        return ns[0].label if hasattr(ns, "__len__") else ns.label
    for inst in odb.rootAssembly.instances.values():
        if nset_name in inst.nodeSets:
            ns = inst.nodeSets[nset_name].nodes[0]
            return ns[0].label if hasattr(ns, "__len__") else ns.label
    raise KeyError(nset_name)

def extract_job(jobname):
    odb_file = jobname + ".odb"
    if not os.path.exists(odb_file):
        print("%-24s | ERRO: ODB nao encontrado" % jobname)
        return None
    Area, F_ref = job_params(jobname)
    odb = openOdb(odb_file, readOnly=True)
    step = odb.steps["PULL"]

    lbl_rp = first_node_label(odb, "RP")
    lbl_ib = first_node_label(odb, "IF_BOT_0")
    lbl_it = first_node_label(odb, "IF_TOP_0")

    u_rp, u_ifb, u_ift = None, None, None
    allwk, allse, alldmd, allvd, lpf_hist = None, None, None, None, None

    for k, hr in step.historyRegions.items():
        ho = hr.historyOutputs
        if "ALLWK" in ho:
            allwk = [x[1] for x in ho["ALLWK"].data]
            allse = [x[1] for x in ho["ALLSE"].data]
            alldmd = [x[1] for x in ho["ALLDMD"].data]
            allvd = [x[1] for x in ho["ALLVD"].data]
            if "LPF" in ho:
                lpf_hist = [x[1] for x in ho["LPF"].data]
        if ("." + str(lbl_rp)) in k or (" " + str(lbl_rp)) in k:
            if "U3" in ho:
                u_rp = ho["U3"].data
        if ("." + str(lbl_ib)) in k or (" " + str(lbl_ib)) in k:
            if "U3" in ho:
                u_ifb = ho["U3"].data
        if ("." + str(lbl_it)) in k or (" " + str(lbl_it)) in k:
            if "U3" in ho:
                u_ift = ho["U3"].data
    odb.close()

    npts = len(u_rp)
    arcs = [x[0] for x in u_rp]
    us = [x[1] for x in u_rp]
    ds = [(u_ift[i][1] - u_ifb[i][1]) for i in range(npts)]
    lpfs = lpf_hist
    tns = [l * F_ref / Area for l in lpfs]

    # Encontrar o primeiro pico fisico (subida inicial com LPF > 0)
    ipk = 0
    for i in range(1, npts):
        if tns[i] > tns[ipk]:
            ipk = i
        if ipk > 5 and tns[i] < 0.85 * tns[ipk]:
            break

    # Apos o pico, parar na rotura completa (Tn <= 0.05 GPa) ou quando ALLDMD deixa de crescer (descarga elastica)
    iend = npts
    stagnant = 0
    for i in range(ipk + 1, npts):
        if tns[i] <= 0.05:
            iend = i + 1
            break
        if alldmd[i] <= alldmd[i - 1] + 1.0e-7 and tns[i] < 0.90 * tns[ipk]:
            stagnant += 1
            if stagnant >= 2:
                iend = i - 1
                break
        else:
            stagnant = 0

    csv_name = jobname + "_curve.csv"
    with open(csv_name, "w") as f:
        f.write("inc,arc,LPF,u_RP_nm,delta_coh_nm,F_nN,Tn_GPa,ALLWK,ALLSE,ALLDMD,ALLVD\n")
        for i in range(iend):
            f.write("%d,%.8f,%.8f,%.8f,%.8f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n" % (
                i, arcs[i], lpfs[i], us[i], ds[i], lpfs[i] * F_ref, tns[i],
                allwk[i], allse[i], alldmd[i], allvd[i]
            ))

    us_c = us[:iend]
    ds_c = ds[:iend]
    tns_c = tns[:iend]
    t_pk = tns_c[ipk]
    u_pk = us_c[ipk]
    u_min_post = min(us_c[ipk:])
    recuo = u_pk - u_min_post
    snap = recuo > 1.0e-3
    G_dmd = alldmd[iend - 1] / Area
    G_vd = allvd[iend - 1] / Area

    print("%-24s | pts=%3d | t_pk=%6.3f GPa @ u_pk=%6.4f nm | recuo=%6.4f nm | d_end=%6.4f nm | Tn_end=%5.3f | G_dmd=%.4f | Snap=%s" % (
        jobname, iend, t_pk, u_pk, recuo, ds_c[-1], tns_c[-1], G_dmd, "SIM" if snap else "NAO"
    ))
    return (jobname, iend, t_pk, u_pk, recuo, ds_c[-1], tns_c[-1], G_dmd, G_vd, "SIM" if snap else "NAO")

if __name__ == "__main__":
    rows = []
    for j in sys.argv[1:]:
        r = extract_job(j)
        if r is not None:
            rows.append(r)
    if rows:
        with open("abaqus_summary.csv", "w") as f:
            f.write("job,pts,t_pk_GPa,u_pk_nm,recuo_nm,d_end_nm,Tn_end_GPa,G_dmd_Jm2,G_vd_Jm2,snap_back\n")
            for r in rows:
                f.write("%s,%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%s\n" % r)

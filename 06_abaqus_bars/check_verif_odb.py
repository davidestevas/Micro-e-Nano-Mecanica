from __future__ import print_function, division
from odbAccess import openOdb

def get_nset(odb, name):
    if name in odb.rootAssembly.nodeSets:
        return odb.rootAssembly.nodeSets[name]
    for inst in odb.rootAssembly.instances.values():
        if name in inst.nodeSets:
            return inst.nodeSets[name]
    raise KeyError("NodeSet not found: " + name)

def get_energy_hist(step):
    for hr in step.historyRegions.values():
        if "ALLWK" in hr.historyOutputs:
            return hr.historyOutputs
    raise KeyError("ALLWK not found in historyRegions: " + str(step.historyRegions.keys()))

for name, expected_t0, expected_Gc in [
    ("verif_coh_bilin", 15.6991, 2.126821),
    ("verif_coh_xn",    15.6991, 2.126821),
]:
    odb = openOdb(name + ".odb", readOnly=True)
    step = odb.steps["PULL"]
    top_set = get_nset(odb, "TOP")
    us, ts = [], []
    for fr in step.frames:
        u3 = fr.fieldOutputs["U"].getSubset(region=top_set).values[0].data[2]
        rf3 = sum(v.data[2] for v in fr.fieldOutputs["RF"].getSubset(region=top_set).values)
        us.append(u3)
        ts.append(rf3)
    G_int = sum(0.5 * (ts[i] + ts[i + 1]) * (us[i + 1] - us[i]) for i in range(len(us) - 1))
    hist = get_energy_hist(step)
    allwk = hist["ALLWK"].data[-1][1]
    alldmd = hist["ALLDMD"].data[-1][1]
    allse = hist["ALLSE"].data[-1][1]
    allvd = hist["ALLVD"].data[-1][1]
    odb.close()
    t_pk = max(ts)
    u_pk = us[ts.index(t_pk)]
    print("=== VERIFICACAO: %s ===" % name)
    print("  Pico T_n  = %.4f GPa (alvo %.4f GPa, erro = %+.3f%%) em delta = %.4f nm" % (t_pk, expected_t0, 100.0 * (t_pk / expected_t0 - 1.0), u_pk))
    print("  Int T du  = %.5f J/m^2 (alvo %.5f J/m^2, erro = %+.3f%%)" % (G_int, expected_Gc, 100.0 * (G_int / expected_Gc - 1.0)))
    print("  Energias  : ALLWK=%.5f | ALLDMD=%.5f | ALLSE=%.5f | ALLVD=%.5f\n" % (allwk, alldmd, allse, allvd))

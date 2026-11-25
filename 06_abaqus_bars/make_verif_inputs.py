from __future__ import print_function, division
import math

C11, C12, C44 = 169.92, 122.63, 76.21
s10 = math.sqrt(10.0)
ax1, ay1, az1 = -1.0 / s10, 0.0, 3.0 / s10
bx1, by1, bz1 =  3.0 / s10, 0.0, 1.0 / s10

# 1) Cubo unitario C3D8 rodado para [310] (verificacao C6)
inp_cube = """*HEADING
Verificacao C3D8 Cu anisotropico rodado para [310] (0 K, unidades nm-nN-GPa)
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 1.0
6, 1.0, 0.0, 1.0
7, 1.0, 1.0, 1.0
8, 0.0, 1.0, 1.0
*ELEMENT, TYPE=C3D8, ELSET=EALL
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=BOT
1, 2, 3, 4
*NSET, NSET=TOP
5, 6, 7, 8
*ORIENTATION, NAME=ORI_G1, DEFINITION=COORDINATES
%.12f, %.12f, %.12f, %.12f, %.12f, %.12f
3, 0.0
*SOLID SECTION, ELSET=EALL, MATERIAL=CU_0K, ORIENTATION=ORI_G1
*MATERIAL, NAME=CU_0K
*ELASTIC, TYPE=ORTHOTROPIC
%.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f, %.4f
%.4f
*STEP, NAME=UNIAX_STRAIN
*STATIC
1.0, 1.0, 1e-5, 1.0
*BOUNDARY
BOT, 1, 3, 0.0
TOP, 1, 2, 0.0
TOP, 3, 3, 1.0e-3
*OUTPUT, FIELD, VARIABLE=PRESELECT
*OUTPUT, HISTORY, VARIABLE=PRESELECT
*END STEP
""" % (ax1, ay1, az1, bx1, by1, bz1, C11, C12, C11, C12, C12, C11, C44, C44, C44)

with open("verif_c3d8_ortho.inp", "w") as f:
    f.write(inp_cube)

# 2) Elemento coesivo COH3D8 unitario (Bilinear 0 K)
t0 = 15.6991
Gc = 2.126821
Kn_bil = 545.70

inp_bilin = """*HEADING
Verificacao COH3D8 Bilinear 0 K (Cu Sigma5(310), unidades nm-nN-GPa)
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 0.0
6, 1.0, 0.0, 0.0
7, 1.0, 1.0, 0.0
8, 0.0, 1.0, 0.0
*ELEMENT, TYPE=COH3D8, ELSET=COH
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=BOT
1, 2, 3, 4
*NSET, NSET=TOP
5, 6, 7, 8
*COHESIVE SECTION, ELSET=COH, MATERIAL=MAT_BILIN, RESPONSE=TRACTION SEPARATION, THICKNESS=SPECIFIED
1.0
*MATERIAL, NAME=MAT_BILIN
*ELASTIC, TYPE=TRACTION
%.4f, %.4f, %.4f
*DAMAGE INITIATION, CRITERION=MAXS
%.6f, %.6f, %.6f
*DAMAGE EVOLUTION, TYPE=ENERGY, SOFTENING=LINEAR
%.6f
*STEP, NAME=PULL, NLGEOM=NO, INC=1000
*STATIC
0.005, 1.0, 1e-6, 0.005
*BOUNDARY
BOT, 1, 3, 0.0
TOP, 1, 2, 0.0
TOP, 3, 3, 0.35
*OUTPUT, FIELD, FREQUENCY=1
*NODE OUTPUT
U, RF
*ELEMENT OUTPUT, ELSET=COH
S, E, SDEG
*OUTPUT, HISTORY, FREQUENCY=1
*ENERGY OUTPUT
ALLWK, ALLSE, ALLDMD, ALLVD
*END STEP
""" % (Kn_bil, Kn_bil, Kn_bil, t0, t0, t0, Gc)

with open("verif_coh_bilin.inp", "w") as f:
    f.write(inp_bilin)

# 3) Elemento coesivo COH3D8 unitario (Xu-Needleman tabular exato 0 K)
dc = Gc / (math.e * t0)
Kn_xn = math.e * t0 / dc
d_init = 0.005 * dc
t_init = Kn_xn * d_init

lines_tab = ["0.00000000, 0.00000000"]
for i in range(1, 121):
    s = (i / 120.0) ** 1.6
    delta = d_init + s * (10.0 * dc - d_init)
    D_val = 1.0 - math.exp(-delta / dc) if i < 120 else 1.0
    u_eff = delta - d_init
    lines_tab.append("%.8f, %.8f" % (D_val, u_eff))
tab_str = "\n".join(lines_tab)

inp_xn = """*HEADING
Verificacao COH3D8 Xu-Needleman Tabular 0 K (Cu Sigma5(310), unidades nm-nN-GPa)
*NODE
1, 0.0, 0.0, 0.0
2, 1.0, 0.0, 0.0
3, 1.0, 1.0, 0.0
4, 0.0, 1.0, 0.0
5, 0.0, 0.0, 0.0
6, 1.0, 0.0, 0.0
7, 1.0, 1.0, 0.0
8, 0.0, 1.0, 0.0
*ELEMENT, TYPE=COH3D8, ELSET=COH
1, 1, 2, 3, 4, 5, 6, 7, 8
*NSET, NSET=BOT
1, 2, 3, 4
*NSET, NSET=TOP
5, 6, 7, 8
*COHESIVE SECTION, ELSET=COH, MATERIAL=MAT_XN, RESPONSE=TRACTION SEPARATION, THICKNESS=SPECIFIED
1.0
*MATERIAL, NAME=MAT_XN
*ELASTIC, TYPE=TRACTION
%.4f, %.4f, %.4f
*DAMAGE INITIATION, CRITERION=MAXS
%.6f, %.6f, %.6f
*DAMAGE EVOLUTION, TYPE=DISPLACEMENT, SOFTENING=TABULAR
%s
*STEP, NAME=PULL, NLGEOM=NO, INC=1000
*STATIC
0.005, 1.0, 1e-6, 0.005
*BOUNDARY
BOT, 1, 3, 0.0
TOP, 1, 2, 0.0
TOP, 3, 3, 0.55
*OUTPUT, FIELD, FREQUENCY=1
*NODE OUTPUT
U, RF
*ELEMENT OUTPUT, ELSET=COH
S, E, SDEG
*OUTPUT, HISTORY, FREQUENCY=1
*ENERGY OUTPUT
ALLWK, ALLSE, ALLDMD, ALLVD
*END STEP
""" % (Kn_xn, Kn_xn, Kn_xn, t_init, t_init, t_init, tab_str)

with open("verif_coh_xn.inp", "w") as f:
    f.write(inp_xn)

print("OK - Gerados: verif_c3d8_ortho.inp, verif_coh_bilin.inp, verif_coh_xn.inp")

from make_bar_inp import write_bar_inp

# 1. Fechar bar0K_L10nm_bilin e bar0K_L40nm_bilin
write_bar_inp("bar0K_L10nm_bilin", 10.0, law="bilin", temp_K=0, riks_inc=600, arc_max=0.012, frac_fillet=0.20)
write_bar_inp("bar0K_L40nm_bilin", 40.0, law="bilin", temp_K=0, riks_inc=800, arc_max=0.012, frac_fillet=0.20)

# 2. Estudo GCI com r = 2.0 CONSTANTE (h = 0.50, 0.25, 0.125 nm) e malha UNIFORME em z (L = 5 nm, 0 K, Xu-Needleman)
for bc in ["pbc", "free"]:
    write_bar_inp("gci_%s_r2_M1" % bc, 5.0, law="xn", temp_K=0, nx=2, ny=3,  h_fine=0.500, z_fine=2.5, h_coarse=0.500, bc=bc, riks_inc=450, arc_max=0.015)
    write_bar_inp("gci_%s_r2_M2" % bc, 5.0, law="xn", temp_K=0, nx=4, ny=6,  h_fine=0.250, z_fine=2.5, h_coarse=0.250, bc=bc, riks_inc=450, arc_max=0.015)
    write_bar_inp("gci_%s_r2_M3" % bc, 5.0, law="xn", temp_K=0, nx=8, ny=12, h_fine=0.125, z_fine=2.5, h_coarse=0.125, bc=bc, riks_inc=450, arc_max=0.015)

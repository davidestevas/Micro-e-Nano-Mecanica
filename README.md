# Multiscale Study of Copper Sigma5(310)[001] Bicrystal Fracture (0 K & 300 K)

Projeto de **Micro e Nano Mecânica** (Mestrado em Engenharia Mecânica, Instituto Superior Técnico — Universidade de Lisboa).

Estudo multiescala das propriedades elásticas e de fratura intergranular de bicristais de cobre (Cu) com fronteira de grão simétrica $\Sigma 5(310)[001]$ a **0 K** e **300 K**, ligando simulações atomísticas (**LAMMPS**, potencial EAM Mishin 2001) a modelos de Elementos Finitos com elementos coesivos (**Abaqus/Standard 2023**, leis Bilinear e Xu–Needleman + método de Riks para captura de *snap-back*).

## Estrutura do Repositório

```text
.
├── 00_potentials/       # Potencial interatómico EAM Cu (Mishin et al., 2001)
├── 01_lammps_bulk/      # Constantes elásticas C11, C12, C44 (0 K e 300 K), NPT e verificação NVE
├── 02_lammps_gb_310/    # Superfície (310) e fronteira Sigma5(310), grid 20x20 e convergência
├── 03_lammps_czm_pull/  # Extração das leis tração-separação (CZM) a 0 K (MS) e 300 K (MD, 1 e 10 m/s)
├── 04_lammps_bars/      # Tração de barras bicristal (L = 2, 5, 10, 20, 40 nm) a 0 K e 300 K
├── 05_python_bridge/    # Scripts Python de análise, calibração CZM (JSON) e geração das figuras finais
├── 06_abaqus_bars/      # Modelos 3D Abaqus 2023 (.inp, .py, curvas .csv, .sta e script PowerShell)
├── 07_report_figs/      # Figuras geradas (.png / .pdf)
└── requirements.txt     # Dependências Python
```

## Pré-requisitos

* **Python 3.9+** (`pip install -r requirements.txt`)
* **LAMMPS** (`lmp_serial` ou `lmp`, com pacote `MANYBODY` para `eam/alloy`)
* **Abaqus 2023** (Abaqus/Standard + `abaqus python`, utilizado em `06_abaqus_bars/`)

## Como Reproduzir o Projeto Passo a Passo

### 1. Propriedades Elásticas do Cu Bulk (`01_lammps_bulk/`)
```bash
cd 01_lammps_bulk
lmp_serial -in in.cu_bulk_0K -log log.cu_bulk_0K
lmp_serial -in in.cu_delta_check_0K -log log.cu_delta_check_0K
lmp_serial -in in.cu_npt_300K -log log.cu_npt_300K
lmp_serial -in in.cu_nve_check -log log.cu_nve_check
lmp_serial -in in.cu_elastic_300K -log log.cu_elastic_300K
lmp_serial -in in.cu_c44_300K -log log.cu_c44_300K
python3 ../05_python_bridge/analyze_npt_300K.py
python3 ../05_python_bridge/check_nve_drift.py
python3 ../05_python_bridge/analyze_elastic_300K.py
```

### 2. Superfície (310) e Fronteira de Grão Sigma5(310) (`02_lammps_gb_310/`)
```bash
cd 02_lammps_gb_310
python3 build_310_structures.py
lmp_serial -in in.cu_surf_310 -log log.cu_surf_310
lmp_serial -in in.cu_gb_grid -log log.cu_gb_grid
python3 ../05_python_bridge/analyze_gb_grid.py
lmp_serial -in in.cu_gb_verify -log log.cu_gb_verify
```

### 3. Extração das Leis Coesivas (CZM) a 0 K e 300 K (`03_lammps_czm_pull/`)
```bash
cd 03_lammps_czm_pull
for H in 5 8 11 15 25; do
  lmp_serial -var Hmob $H -in in.cu_pull_0K -log log.pull0K_Hmob$(printf "%02d" $H)
done
python3 ../05_python_bridge/analyze_czm_0K.py

lmp_serial -in in.cu_gb_eq_300K_v2 -log log.gb_eq_300K_v2
lmp_serial -var vel 0.10 -var nsteps 100000 -var ave_n 50  -var tag v10ms -in in.cu_pull_300K -log log.pull300K_v10ms
lmp_serial -var vel 0.01 -var nsteps 1000000 -var ave_n 500 -var tag v1ms  -in in.cu_pull_300K -log log.pull300K_v1ms
python3 ../05_python_bridge/analyze_czm_300K.py
```

### 4. Tração das Barras Bicristal em LAMMPS (`04_lammps_bars/`)
```bash
cd 04_lammps_bars
python3 build_bars.py

# 0 K (MS quase-estático)
lmp_serial -var Lnm 2  -var du 0.02 -var nsteps 400 -var dskip 50 -in in.bar_0K_pull -log log.bar0K_L2nm
lmp_serial -var Lnm 5  -var du 0.02 -var nsteps 550 -var dskip 50 -in in.bar_0K_pull -log log.bar0K_L5nm
lmp_serial -var Lnm 10 -var du 0.02 -var nsteps 750 -var dskip 75 -in in.bar_0K_pull -log log.bar0K_L10nm
lmp_serial -var Lnm 20 -var du 0.04 -var nsteps 600 -var dskip 60 -in in.bar_0K_pull -log log.bar0K_L20nm
lmp_serial -var Lnm 40 -var du 0.05 -var nsteps 840 -var dskip 84 -in in.bar_0K_pull -log log.bar0K_L40nm
python3 extract_last_frames.py
lmp_serial -var Lnm 20 -var u0 24.0 -var dE0 22.0653 -var du 0.05 -var nsteps 120 -var dskip 20 -in in.bar_0K_cont -log log.bar0K_L20nm_cont
lmp_serial -var Lnm 40 -var u0 42.0 -var dE0 34.9780 -var du 0.10 -var nsteps 150 -var dskip 25 -in in.bar_0K_cont -log log.bar0K_L40nm_cont
python3 ../05_python_bridge/analyze_bars_0K.py

# 300 K (MD NVT a 10 m/s e 1 m/s + verificação de seed)
lmp_serial -var Lnm 5  -var vel 0.10 -var nsteps 120000  -var ave_n 50  -var seed 492845 -var tag L5nm_v10ms    -in in.bar_300K_pull -log log.bar300K_L5nm_v10ms
lmp_serial -var Lnm 5  -var vel 0.10 -var nsteps 120000  -var ave_n 50  -var seed 918273 -var tag L5nm_v10ms_s2 -in in.bar_300K_pull -log log.bar300K_L5nm_v10ms_s2
lmp_serial -var Lnm 5  -var vel 0.01 -var nsteps 1200000 -var ave_n 500 -var seed 492845 -var tag L5nm_v1ms     -in in.bar_300K_pull -log log.bar300K_L5nm_v1ms
lmp_serial -var Lnm 10 -var vel 0.10 -var nsteps 180000  -var ave_n 50  -var seed 492845 -var tag L10nm_v10ms   -in in.bar_300K_pull -log log.bar300K_L10nm_v10ms
lmp_serial -var Lnm 20 -var vel 0.10 -var nsteps 280000  -var ave_n 100 -var seed 492845 -var tag L20nm_v10ms   -in in.bar_300K_pull -log log.bar300K_L20nm_v10ms
lmp_serial -var Lnm 40 -var vel 0.10 -var nsteps 480000  -var ave_n 150 -var seed 492845 -var tag L40nm_v10ms   -in in.bar_300K_pull -log log.bar300K_L40nm_v10ms
python3 extract_cont_300K.py
lmp_serial -var Lnm 20 -var nsteps 105000 -var ave_n 100 -in in.bar_300K_cont -log log.bar300K_L20nm_cont
lmp_serial -var Lnm 40 -var nsteps 150000 -var ave_n 150 -in in.bar_300K_cont -log log.bar300K_L40nm_cont
python3 ../05_python_bridge/analyze_bars_300K.py
```

### 5. Elementos Finitos 3D no Abaqus 2023 (`06_abaqus_bars/`, Windows PowerShell)
Sistema de unidades consistentes: **nm, nN, GPa** ($1\;\mathrm{J/m^2} = 1\;\mathrm{GPa\cdot nm}$).
```powershell
cd 06_abaqus_bars
.\run_all_abaqus.ps1
```

### 6. Figuras Comparativas Finais LAMMPS vs. Abaqus e Métricas GCI
```bash
python3 05_python_bridge/make_final_figures.py
```

## Principais Resultados

| Propriedade / Métrica | 0 K | 300 K (10 m/s) | 300 K (1 m/s) |
| :--- | :---: | :---: | :---: |
| Módulo axial efetivo $M_{\text{free-shear}}^{[310]}$ (GPa) | 177.73 | 171.19 | 171.19 |
| Tração coesiva de pico $t_0$ (GPa) | 15.70 | 15.91 | 15.03 |
| Energia de separação $G_c$ (J/m$^2$) | 2.1268 | 3.2449 | 3.0626 |
| Limiar de *snap-back* Xu–Needleman $L_{\text{snap}}$ (nm) | 1.53 | 2.20 | 2.32 |
| Limiar de *snap-back* Bilinear $L_{\text{snap}}$ (nm) | 2.74 | 3.89 | 4.21 |
| GCI em 3 malhas ($L = 5\;\text{nm}$, 0 K, Xu–Needleman) | 0.076% (PBC) / 0.009% (livre) | — | — |

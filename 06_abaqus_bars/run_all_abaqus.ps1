# Script de reproducao completa da fase Abaqus 2023 (executar dentro de 06_abaqus_bars/)
$ErrorActionPreference = "Continue"

function Run-Abq($job) {
    for ($try = 1; $try -le 3; $try++) {
        Write-Host "`n>>> [$job] tentativa$try..." -ForegroundColor Cyan
        abaqus job=$job interactive ask_delete=OFF | Out-Host
        if (Test-Path "$job.sta") {
            if ((Get-Content "$job.sta" -Tail 5 | Out-String) -match "COMPLETED SUCCESSFULLY") { return }
        }
        Start-Sleep -Seconds 5
    }
    Write-Host "AVISO: $job nao completou apos 3 tentativas." -ForegroundColor Yellow
}

# 1. Verificacao de codigo (1 elemento C3D8 [310] + 1 elemento COH3D8 Bilinear e Xu-Needleman)
abaqus python make_verif_inputs.py
foreach ($v in @("verif_c3d8_ortho", "verif_coh_bilin", "verif_coh_xn")) { Run-Abq $v }
abaqus python check_verif_odb.py

# 2. Barras bicristal 3D (0 K, 300 K e estudo de malha GCI)
$allJobs = @(
    "bar0K_L1nm_bilin", "bar0K_L1nm_xn",
    "bar0K_L2nm_bilin", "bar0K_L2nm_xn",
    "bar0K_L3nm_bilin", "bar0K_L3nm_xn",
    "bar0K_L5nm_bilin", "bar0K_L5nm_xn",
    "bar0K_L10nm_bilin", "bar0K_L10nm_xn",
    "bar0K_L20nm_bilin", "bar0K_L20nm_xn",
    "bar0K_L40nm_bilin", "bar0K_L40nm_xn",
    "bar300K_L2nm_bilin", "bar300K_L2nm_xn",
    "bar300K_L3nm_bilin", "bar300K_L3nm_xn",
    "bar300K_L5nm_bilin", "bar300K_L5nm_xn",
    "bar300K_L10nm_xn", "bar300K_L20nm_xn", "bar300K_L40nm_xn",
    "bar300K1ms_L5nm_bilin", "bar300K1ms_L5nm_xn",
    "gci_pbc_M1", "gci_pbc_M3",
    "gci_free_M1", "gci_free_M2", "gci_free_M3"
)
foreach ($j in $allJobs) { Run-Abq$j }

# 3. Extracao das curvas CSV e tabela resumo
abaqus python extract_bar_odb.py @allJobs

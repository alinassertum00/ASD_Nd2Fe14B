# Portable command examples

All commands are run from the repository root unless stated otherwise.

## 1. Install Python dependencies
```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

## 2. Tell the repository where VAMPIRE is
Either put `vampire-serial` / `vampire-parallel` on PATH, pass `--vampire`, or set:
```bash
export VAMPIRE_BIN=/absolute/path/to/vampire-serial
```
PowerShell:
```powershell
$env:VAMPIRE_BIN="C:\path\to\vampire-serial.exe"
```

## 3. Generate a 50 K CMC sweep
```bash
python simulations/constrained_monte_carlo/generate_cmc_sweep.py \
  --temperature 50 --step-angle 3 --output work/CMC_50K
```

## 4. Run the generated CMC sweep
```bash
python scripts/run_cmc_sweep.py work/CMC_50K --jobs 8
```

## 5. Run one domain-wall case
```bash
python scripts/run_vampire.py simulations/domain_wall
```

## 6. Process a domain-wall output
```bash
python analysis/analyze_dw_profile.py \
  --dw-file work/DW_400K/dw-600000.txt \
  --input-file work/DW_400K/input \
  --output-dir results/400K
```

## 7. Fit a Bloch wall
```bash
python analysis/fit_bloch_wall.py results/400K/dw_profile_step_600000.csv \
  --output-prefix results/400K/bloch
```

## 8. Fit exchange stiffness
```bash
python analysis/fit_exchange_stiffness.py results/50K/dw_profile_step_600000.csv \
  --temperature 50 --k-file config/anisotropy_constants_final.csv \
  --output-prefix results/50K/Ae
```

## 9. Reconstruct exchange shells
```bash
python analysis/shell_pair_excel_export.py model/Nd2Fe14B_scaled_175.ucf \
  --output results/Nd2Fe14B_shell_pairs.xlsx
```

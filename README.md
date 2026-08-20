# ASD_Nd2Fe14B

Portable simulation and post-processing workflow accompanying the thesis
**Finite-Temperature Atomistic Spin Modelling of Intrinsic Magnetic Properties in Nd2Fe14B**.

The repository separates the physical model, VAMPIRE simulation templates, numerical analysis,
and legacy/debugging material. No maintained script assumes a username, home directory, VM path,
or fixed VAMPIRE installation location.

## Repository tree

```text
ASD_Nd2Fe14B/
├── model/                         # Canonical MAT/UCF model files
├── simulations/
│   ├── curie_temperature/         # M(T) / Curie-temperature input
│   ├── constrained_monte_carlo/   # CMC sweep generator + runner
│   └── domain_wall/               # LLG domain-wall input/template
├── analysis/
│   ├── fit_anisotropy_direct_torque.py
│   ├── analyze_dw_profile.py
│   ├── fit_bloch_wall.py
│   ├── fit_exchange_stiffness.py
│   ├── shell_pair_excel_export.py
│   └── validate_model_files.py
├── scripts/
│   ├── run_vampire.py             # Portable single-case runner
│   └── run_cmc_sweep.py           # Portable parallel CMC runner
├── config/
│   └── anisotropy_constants_final.csv
├── examples/
│   └── portable_commands.md
├── results/                        # Generated results; ignored by git
├── archive/                        # Original one-off/debug scripts
├── docs/
├── requirements.txt
└── .gitignore
```

## Requirements

1. Python 3.10+ recommended.
2. A working VAMPIRE installation. VAMPIRE itself is **not bundled** in this repository.
3. Python packages:

```bash
python -m venv .venv
source .venv/bin/activate              # Linux/macOS
# .venv\Scripts\Activate.ps1          # Windows PowerShell
python -m pip install -r requirements.txt
```

## Connecting VAMPIRE on any machine

The maintained Python runners find the executable in this order:

1. explicit `--vampire /path/to/vampire-serial`;
2. environment variable `VAMPIRE_BIN`;
3. `vampire-serial` available on `PATH`.

Linux/macOS example:
```bash
export VAMPIRE_BIN=/opt/vampire/vampire-serial
```

Windows/PowerShell example (for a compatible executable/environment):
```powershell
$env:VAMPIRE_BIN="C:\path\to\vampire-serial.exe"
```

No source file needs to be edited to change machines.

## Core workflows

### Curie temperature / magnetisation
```bash
python scripts/run_vampire.py simulations/curie_temperature --vampire /path/to/vampire-parallel
```

### Constrained Monte Carlo anisotropy
Generate a sweep using model files from `model/` automatically:
```bash
python simulations/constrained_monte_carlo/generate_cmc_sweep.py \
  --temperature 50 --step-angle 3 --output work/CMC_50K
```
Run it:
```bash
python scripts/run_cmc_sweep.py work/CMC_50K --vampire /path/to/vampire-serial --jobs 8
```
Fit K1-K3 after collecting the CMC output:
```bash
python analysis/fit_anisotropy_direct_torque.py collected_50K.xlsx \
  --auto-filter --volume-m3 1.2e-24 --output-prefix results/50K_K
```

### Domain wall
Run a prepared domain-wall directory:
```bash
python scripts/run_vampire.py work/DW_400K --vampire /path/to/vampire-serial
```
Build the 1-D plane-averaged profile:
```bash
python analysis/analyze_dw_profile.py \
  --dw-file work/DW_400K/dw-600000.txt \
  --input-file work/DW_400K/input \
  --output-dir results/400K
```
Fit Bloch width:
```bash
python analysis/fit_bloch_wall.py results/400K/dw_profile_step_600000.csv \
  --output-prefix results/400K/bloch
```
Fit exchange stiffness using the full K1-K2-K3 wall relation:
```bash
python analysis/fit_exchange_stiffness.py results/50K/dw_profile_step_600000.csv \
  --temperature 50 --k-file config/anisotropy_constants_final.csv \
  --output-prefix results/50K/Ae
```

### Exchange-shell documentation
```bash
python analysis/shell_pair_excel_export.py model/Nd2Fe14B_scaled_175.ucf \
  --output results/Nd2Fe14B_shell_pairs.xlsx
```

## Important scientific conventions

- `fit_bloch_wall.py` is intended for Bloch-like walls. Do not assign the broad low-temperature
  theta0 -> -theta0 metastable wall a conventional Bloch width with this model.
- `fit_exchange_stiffness.py` retains K1, K2 and K3 and uses the general wall integral.
- Processed `x_minus_x0_nm` from `analyze_dw_profile.py` is already in physical nanometres because
  its spacing is obtained from `Lx / Nbins`; therefore the exchange-stiffness fit defaults to
  `--x-scale 1.0`.
- The torque-density conversion volume is an explicit CLI parameter. The thesis-calibration default
  in `fit_anisotropy_direct_torque.py` is 1.2e-24 m^3; users applying the script to another geometry
  should provide their appropriate value.

## Legacy files

`archive/original_uploaded_scripts/` contains historical debugging and machine-specific scripts for
provenance. They are intentionally not part of the maintained workflow.

## More examples
See `examples/portable_commands.md` and `docs/CODE_REVIEW.md`.

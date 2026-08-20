# Code review and corrections

## Uploaded Bloch-wall script
**Status:** core fitting equations are correct; packaging was not repository-ready.

Retained:
- `Mz = C - A*tanh((x-x0)/delta)`
- `Mxy = C + A/cosh((x-x0)/delta)`
- `Delta_DW = pi*|delta|`

Corrected:
- removed `cat > ... EOF` shell wrapper from the Python file;
- removed hard-coded profile name;
- added CSV/XLSX input support and column validation;
- exports all fitted parameters and R2 values with explicit names.

## Uploaded K-fit scripts
**Status:** torque model is useful; the files should not be used unchanged.

Problems found:
- both use a hard-coded volume of `1.1e-24 m^3`;
- hard-coded data locations make the result difficult to reproduce on another machine;
- `another_K_fit.txt` mixes a 0 T dataset and a collective 5 T comparison dataset and is
  therefore a comparison/diagnostic script rather than the production fit;
- temperature-smoothness and K3 penalties are hard-coded. Regularisation can be useful for
  the ill-conditioned sixth-order basis, but it must be explicit rather than hidden.

The maintained `fit_anisotropy_direct_torque.py` therefore fits one dataset transparently,
uses the final project volume of `1.2e-24 m^3`, exposes the fitting window and optional K3
penalty, and prints the matrix condition number.

## Uploaded CMC scripts
**Status:** valid diagnostic helper, not a general production script.

The original generator creates only six 10 K cases between 65.172 and 80.690 degrees and uses
50,000 equilibration plus 50,000 averaging steps. The runner launches all angle folders at once.
The maintained generator accepts temperature/angle/step settings as arguments and the runner
limits the number of simultaneous jobs.

## Material file
**Status:** physical site parameters are useful, but the uploaded file ends with exploratory
Case-D domain-wall initialization settings. Those initialization lines are not intrinsic material
parameters and must not live in the canonical MAT used for Curie-temperature or CMC runs.

The project archive separates:
- `Nd2Fe14B_base.mat` — physical material model;
- `Nd2Fe14B_DW_180deg.mat` — standard +z/-z domain-wall initialization;
- `Nd2Fe14B_DW_caseA_50K.mat` — targeted Case-A low-temperature initialization.

## Exchange stiffness
The new `fit_exchange_stiffness.py` uses the full anisotropy landscape and does not silently
replace it by the high-temperature approximation `(Delta/pi)^2*K1`.

### Important coordinate correction
The processed 50 K summary reports `Lx_nm=40` and `dx_nm=1.0`. Therefore a column already
called `x_minus_x0_nm` is a physical nm coordinate. Applying an additional factor of 0.8776004
would double-count a length conversion. The maintained code defaults to `x_scale=1.0` and only
allows another scale when the input coordinate is explicitly an integer unit-cell/bin index.

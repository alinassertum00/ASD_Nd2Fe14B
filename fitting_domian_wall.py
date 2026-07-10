cat > fit_bloch_wall_corrected.py << 'EOF'
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# CHANGE ONLY THIS LINE
csv_file = "dw_profile_step_600000.csv"

df = pd.read_csv(csv_file)

x   = df["x_minus_x0_nm"].values
Mz  = df["Mz"].values
Mxy = df["Mxy"].values

# Correct Bloch-like models with background offsets
def mz_model(x, C, A, delta, x0):
    return C - A*np.tanh((x-x0)/delta)

def mxy_model(x, C, A, delta, x0):
    return C + A/np.cosh((x-x0)/delta)

def r2(y, yf):
    return 1.0 - np.sum((y-yf)**2)/np.sum((y-y.mean())**2)

# Fit Mz
popt_mz, _ = curve_fit(
    mz_model, x, Mz,
    p0=[0.0, (Mz.max()-Mz.min())/2, 1.0, 0.0],
    maxfev=100000
)

# Fit Mxy
popt_mxy, _ = curve_fit(
    mxy_model, x, Mxy,
    p0=[Mxy.min(), Mxy.max()-Mxy.min(), abs(popt_mz[2]), 0.0],
    maxfev=100000
)

Mz_fit  = mz_model(x, *popt_mz)
Mxy_fit = mxy_model(x, *popt_mxy)

delta_mz  = abs(popt_mz[2])
delta_mxy = abs(popt_mxy[2])

dw_mz  = np.pi * delta_mz
dw_mxy = np.pi * delta_mxy
dw_avg = 0.5 * (dw_mz + dw_mxy)

out = pd.DataFrame({
    "x_minus_x0_nm": x,
    "Mz_data": Mz,
    "Mz_fit": Mz_fit,
    "Mxy_data": Mxy,
    "Mxy_fit": Mxy_fit
})
out.to_csv("bloch_fit_corrected_excel_data.csv", index=False)

summary = pd.DataFrame([{
    "csv_file": csv_file,
    "delta_Mz_nm": delta_mz,
    "DW_Mz_pi_delta_nm": dw_mz,
    "R2_Mz": r2(Mz, Mz_fit),
    "delta_Mxy_nm": delta_mxy,
    "DW_Mxy_pi_delta_nm": dw_mxy,
    "R2_Mxy": r2(Mxy, Mxy_fit),
    "Final_average_DW_width_nm": dw_avg
}])
summary.to_csv("bloch_fit_corrected_summary.csv", index=False)

print(summary.T)
print()
print("Saved:")
print("bloch_fit_corrected_excel_data.csv")
print("bloch_fit_corrected_summary.csv")
EOF

python3 fit_bloch_wall_corrected.py

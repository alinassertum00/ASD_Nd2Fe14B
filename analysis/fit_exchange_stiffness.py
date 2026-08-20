#!/usr/bin/env python3
"""Fit the exchange stiffness Ae(T) from a Bloch-like domain-wall profile.

The fit implements the general continuum relation used in the thesis:

    F(theta,T) = K1 sin^2(theta) + K2 sin^4(theta) + K3 sin^6(theta)

    x(theta)-x0 = sqrt(Ae) * integral[dTheta / sqrt(F(Theta)-F(theta0))]

where theta0 is the equilibrium domain angle obtained from the minimum of F.
The integral is referenced to pi/2 to avoid evaluating the endpoint singularity.
Because the processed domain-wall profile is already recentered, the fit is
performed through the origin in x versus the numerical integral I(theta).

Important
---------
* Use the narrow Bloch-like theta0 -> pi-theta0 wall below spin reorientation.
* Do NOT use the broad theta0 -> -theta0 metastable wall with this 1-D polar-angle model.
* K1, K2 and K3 are all retained; K3 is not silently set to zero.
* x_scale converts the profile coordinate to physical distance. In the thesis
  workflow, x_minus_x0_nm is already a physical nanometre coordinate when it was built
  from Lx/Nbins. Therefore the default scale is 1.0. Use a different scale only if
  your input coordinate is explicitly a lattice-cell index rather than nanometres.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import interp1d
from scipy.optimize import minimize_scalar


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported profile format: {path.suffix}")


def free_energy(theta, k1, k2, k3):
    s2 = np.sin(theta) ** 2
    return (k1 * s2 + k2 * s2**2 + k3 * s2**3) * 1e6  # J/m^3


def equilibrium_angle(k1, k2, k3):
    result = minimize_scalar(
        lambda t: free_energy(t, k1, k2, k3),
        bounds=(0.0, 0.5 * np.pi),
        method="bounded",
        options={"xatol": 1e-14},
    )
    return float(result.x)


def build_integral_interpolator(theta0, k1, k2, k3, n_grid=250000):
    # Avoid the exact equilibrium endpoints where the integrand diverges.
    eps = max(1e-7, 1e-6 * max(1.0, theta0))
    low = theta0 + eps
    high = np.pi - theta0 - eps
    if low >= high:
        raise ValueError("Invalid integration interval; check K values.")

    theta_grid = np.linspace(low, high, n_grid)
    f0 = free_energy(theta0, k1, k2, k3)
    delta_f = free_energy(theta_grid, k1, k2, k3) - f0

    # Numerical roundoff can produce tiny negative values near the minimum.
    positive = delta_f[delta_f > 0]
    if len(positive) == 0:
        raise ValueError("Anisotropy landscape does not provide a valid wall integral.")
    floor = max(np.min(positive) * 1e-6, 1e-20)
    integrand = 1.0 / np.sqrt(np.maximum(delta_f, floor))

    cumulative = cumulative_trapezoid(integrand, theta_grid, initial=0.0)
    reference = np.interp(0.5 * np.pi, theta_grid, cumulative)
    integral = cumulative - reference

    return interp1d(
        theta_grid,
        integral,
        kind="linear",
        bounds_error=False,
        fill_value=np.nan,
    )


def load_k_from_csv(path: Path, temperature: float):
    table = pd.read_csv(path)
    required = {"Temperature_K", "K1_MJ_m3", "K2_MJ_m3", "K3_MJ_m3"}
    if not required.issubset(table.columns):
        raise KeyError(f"K table must contain columns: {sorted(required)}")

    distances = np.abs(table["Temperature_K"].to_numpy(dtype=float) - temperature)
    idx = int(np.argmin(distances))
    if distances[idx] > 1e-9:
        raise ValueError(
            f"No exact K row for T={temperature:g} K in {path}. "
            "Provide --k1 --k2 --k3 explicitly rather than interpolating silently."
        )
    row = table.iloc[idx]
    return float(row.K1_MJ_m3), float(row.K2_MJ_m3), float(row.K3_MJ_m3)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="Processed wall profile CSV/XLSX")
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--k-file", type=Path, help="CSV table of K1,K2,K3 versus temperature")
    parser.add_argument("--k1", type=float, help="K1 in MJ/m^3")
    parser.add_argument("--k2", type=float, help="K2 in MJ/m^3")
    parser.add_argument("--k3", type=float, help="K3 in MJ/m^3")
    parser.add_argument("--x-column", default="x_minus_x0_nm")
    parser.add_argument("--mz-column", default="Mz")
    parser.add_argument(
        "--x-scale", type=float, default=1.0,
        help="Physical nm per unit in x column (default: 1.0 for a column already in nm)",
    )
    parser.add_argument(
        "--fit-half-width-nm", type=float, default=2.0,
        help="Use only |x| <= this physical half-width for the stable central fit",
    )
    parser.add_argument("--output-prefix", default="exchange_stiffness_fit")
    args = parser.parse_args()

    if args.k_file is not None:
        k1, k2, k3 = load_k_from_csv(args.k_file, args.temperature)
    else:
        if None in (args.k1, args.k2, args.k3):
            parser.error("Use --k-file or supply --k1, --k2 and --k3.")
        k1, k2, k3 = args.k1, args.k2, args.k3

    df = read_table(args.profile)
    for column in (args.x_column, args.mz_column):
        if column not in df.columns:
            raise KeyError(f"Missing required column: {column}")

    x_file = df[args.x_column].to_numpy(dtype=float)
    mz = df[args.mz_column].to_numpy(dtype=float)
    finite = np.isfinite(x_file) & np.isfinite(mz)
    x_file, mz = x_file[finite], mz[finite]

    # Physical, centered coordinate in metres.
    x_nm = x_file * args.x_scale
    x_m = x_nm * 1e-9

    # For the normalized wall profiles used in the thesis, Mz = cos(theta).
    if np.nanmax(np.abs(mz)) > 1.05:
        raise ValueError(
            "Mz does not look normalized. Normalize by Ms(T) before this fit, "
            "or provide a processed normalized profile."
        )
    theta = np.arccos(np.clip(mz, -1.0, 1.0))

    theta0 = equilibrium_angle(k1, k2, k3)
    integral_of_theta = build_integral_interpolator(theta0, k1, k2, k3)
    i_values = integral_of_theta(theta)

    central = (
        np.isfinite(i_values)
        & (np.abs(x_nm) <= args.fit_half_width_nm)
    )
    if np.count_nonzero(central) < 3:
        raise ValueError(
            "Fewer than three valid central points. Increase --fit-half-width-nm "
            "or verify the profile and x scaling."
        )

    # x = sqrt(Ae) I(theta). The profile has already been shifted by x0,
    # therefore the physically consistent regression is through the origin.
    i_fit = i_values[central]
    x_fit = x_m[central]
    slope = float(np.dot(i_fit, x_fit) / np.dot(i_fit, i_fit))
    ae_j_m = slope**2
    ae_pj_m = ae_j_m * 1e12

    x_model_m = slope * i_values
    residual_nm = (x_m - x_model_m) * 1e9
    rmse_nm = float(np.sqrt(np.nanmean(residual_nm[central] ** 2)))

    prefix = Path(args.output_prefix)
    detail_path = prefix.with_name(prefix.name + "_profile.csv")
    summary_path = prefix.with_name(prefix.name + "_summary.csv")

    pd.DataFrame({
        "x_file": x_file,
        "x_physical_nm": x_nm,
        "Mz": mz,
        "theta_rad": theta,
        "I_theta_sqrt_m3_over_J": i_values,
        "x_model_nm": x_model_m * 1e9,
        "residual_nm": residual_nm,
        "used_in_central_fit": central,
    }).to_csv(detail_path, index=False)

    pd.DataFrame([{
        "Temperature_K": args.temperature,
        "K1_MJ_m3": k1,
        "K2_MJ_m3": k2,
        "K3_MJ_m3": k3,
        "theta0_deg": np.degrees(theta0),
        "x_scale_nm_per_input_unit": args.x_scale,
        "fit_half_width_nm": args.fit_half_width_nm,
        "n_fit_points": int(np.count_nonzero(central)),
        "Ae_J_m": ae_j_m,
        "Ae_pJ_m": ae_pj_m,
        "central_fit_RMSE_nm": rmse_nm,
    }]).to_csv(summary_path, index=False)

    print(f"T = {args.temperature:g} K")
    print(f"K1,K2,K3 = {k1:g}, {k2:g}, {k3:g} MJ/m^3")
    print(f"theta0 = {np.degrees(theta0):.4f} deg")
    print(f"fit points = {np.count_nonzero(central)}")
    print(f"Ae = {ae_pj_m:.6f} pJ/m")
    print(f"central-fit RMSE = {rmse_nm:.6f} nm")
    print(f"Saved: {detail_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()

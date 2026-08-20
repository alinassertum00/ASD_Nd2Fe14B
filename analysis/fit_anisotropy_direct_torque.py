#!/usr/bin/env python3
"""Fit K1, K2 and K3 directly to constrained-Monte-Carlo torque data.

Free-energy density:
    F(theta) = K1 sin^2(theta) + K2 sin^4(theta) + K3 sin^6(theta)

Restoring torque density:
    Ty(theta) = -dF/dtheta
              = -2 K1 sin(theta)cos(theta)
                -4 K2 sin^3(theta)cos(theta)
                -6 K3 sin^5(theta)cos(theta)

The input torque produced by VAMPIRE is a total torque in joules. It is converted
into a torque density using the calibration volume adopted in the final thesis
workflow, V = 1.1e-24 m^3. The volume is exposed as a command-line option so the
assumption is never hidden.

The script can automatically exclude the unstable high-angle portion of a CMC
sweep using the magnetisation-length drop. Because the sin^2/sin^4/sin^6 basis
is strongly correlated, the reported condition number should always be checked.
No literature K values are used as fitting targets.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported input format: {path.suffix}")


def find_column(df, requested, candidates):
    if requested:
        if requested not in df.columns:
            raise KeyError(f"Column '{requested}' not found. Available: {list(df.columns)}")
        return requested
    for name in candidates:
        if name in df.columns:
            return name
    return None


def design_matrix(theta):
    s, c = np.sin(theta), np.cos(theta)
    return np.column_stack((
        -2.0 * s * c,
        -4.0 * s**3 * c,
        -6.0 * s**5 * c,
    ))


def energy(theta, k):
    k1, k2, k3 = k
    s2 = np.sin(theta) ** 2
    return k1 * s2 + k2 * s2**2 + k3 * s2**3


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data", type=Path, help="Collected CMC data CSV/XLSX")
    p.add_argument("--angle-column", default=None)
    p.add_argument("--mag-column", default=None)
    p.add_argument("--torque-column", default=None)
    p.add_argument("--volume-m3", type=float, default=1.1e-24)
    p.add_argument("--torque-is-density", action="store_true",
                   help="Input torque column is already in MJ/m^3")
    p.add_argument("--max-angle-deg", type=float, default=None,
                   help="Explicit maximum angle to include")
    p.add_argument("--auto-filter", action="store_true",
                   help="Stop before first post-30deg M drop below threshold")
    p.add_argument("--m-threshold", type=float, default=0.95,
                   help="Fraction of maximum M used by --auto-filter")
    p.add_argument("--minimum-retained-angle", type=float, default=57.0)
    p.add_argument("--ridge-k3", type=float, default=0.0,
                   help="Optional Tikhonov penalty on K3; default 0 (no regularisation)")
    p.add_argument("--output-prefix", default="anisotropy_fit")
    args = p.parse_args()

    df = read_table(args.data)

    angle_col = find_column(df, args.angle_column,
                            ["phi_deg", "Angle_deg", "angle_deg", "constraint_phi", "Phi_deg"])
    mag_col = find_column(df, args.mag_column,
                          ["M", "mean_magnetisation_length", "Magnetisation", "Output_col_04"])
    torque_col = find_column(df, args.torque_column,
                             ["Ty_J", "Ty", "mean_total_torque_y", "Torque_y_J", "Output_col_06"])

    # Compatibility with the spreadsheets used in the original working scripts:
    # column D = angle, E = |M|, G = Ty. This fallback is explicit and reported.
    used_positional_fallback = False
    if angle_col is None or mag_col is None or torque_col is None:
        if df.shape[1] >= 7:
            angle = pd.to_numeric(df.iloc[:, 3], errors="coerce").to_numpy()
            magnetisation = pd.to_numeric(df.iloc[:, 4], errors="coerce").to_numpy()
            torque_raw = pd.to_numeric(df.iloc[:, 6], errors="coerce").to_numpy()
            used_positional_fallback = True
        else:
            raise KeyError(
                "Could not detect angle/M/torque columns. Specify them explicitly with "
                "--angle-column, --mag-column and --torque-column."
            )
    else:
        angle = pd.to_numeric(df[angle_col], errors="coerce").to_numpy()
        magnetisation = pd.to_numeric(df[mag_col], errors="coerce").to_numpy()
        torque_raw = pd.to_numeric(df[torque_col], errors="coerce").to_numpy()

    valid = np.isfinite(angle) & np.isfinite(magnetisation) & np.isfinite(torque_raw)
    angle, magnetisation, torque_raw = angle[valid], magnetisation[valid], torque_raw[valid]
    order = np.argsort(angle)
    angle, magnetisation, torque_raw = angle[order], magnetisation[order], torque_raw[order]

    if args.torque_is_density:
        torque_density = torque_raw.astype(float)
    else:
        torque_density = torque_raw / args.volume_m3 * 1e-6  # J/m^3 -> MJ/m^3

    include = np.ones_like(angle, dtype=bool)
    if args.max_angle_deg is not None:
        include &= angle <= args.max_angle_deg

    if args.auto_filter:
        mmax = float(np.max(magnetisation))
        bad = np.where((angle >= 30.0) & (magnetisation < args.m_threshold * mmax))[0]
        cutoff = len(angle) - 1 if len(bad) == 0 else max(0, int(bad[0] - 1))
        eligible = np.where(angle <= args.minimum_retained_angle)[0]
        if len(eligible):
            cutoff = max(cutoff, int(eligible[-1]))
        include &= np.arange(len(angle)) <= cutoff

    theta = np.deg2rad(angle[include])
    y = torque_density[include]
    xmat = design_matrix(theta)

    # Optional K3 stabilization only when explicitly requested.
    if args.ridge_k3 > 0:
        penalty = np.array([[0.0, 0.0, np.sqrt(args.ridge_k3)]])
        xsolve = np.vstack([xmat, penalty])
        ysolve = np.concatenate([y, [0.0]])
    else:
        xsolve, ysolve = xmat, y

    k, _, _, _ = np.linalg.lstsq(xsolve, ysolve, rcond=None)
    fitted_all = design_matrix(np.deg2rad(angle)) @ k
    energy_all = energy(np.deg2rad(angle), k)
    rmse = float(np.sqrt(np.mean((y - xmat @ k) ** 2)))
    condition = float(np.linalg.cond(xmat))

    prefix = Path(args.output_prefix)
    detail_path = prefix.with_name(prefix.name + "_all_data.csv")
    summary_path = prefix.with_name(prefix.name + "_summary.csv")

    pd.DataFrame({
        "angle_deg": angle,
        "magnetisation_length": magnetisation,
        "torque_input": torque_raw,
        "torque_density_MJ_m3": torque_density,
        "torque_fit_MJ_m3": fitted_all,
        "free_energy_fit_MJ_m3": energy_all,
        "included_in_fit": include,
    }).to_csv(detail_path, index=False)

    pd.DataFrame([{
        "input_file": str(args.data),
        "volume_m3": args.volume_m3,
        "K1_MJ_m3": k[0],
        "K2_MJ_m3": k[1],
        "K3_MJ_m3": k[2],
        "fit_RMSE_MJ_m3": rmse,
        "design_matrix_condition_number": condition,
        "n_points": int(np.count_nonzero(include)),
        "max_included_angle_deg": float(np.max(angle[include])),
        "ridge_k3": args.ridge_k3,
        "used_positional_D_E_G_fallback": used_positional_fallback,
    }]).to_csv(summary_path, index=False)

    print(f"K1 = {k[0]:.6f} MJ/m^3")
    print(f"K2 = {k[1]:.6f} MJ/m^3")
    print(f"K3 = {k[2]:.6f} MJ/m^3")
    print(f"RMSE = {rmse:.6f} MJ/m^3")
    print(f"Condition number = {condition:.3e}")
    if condition > 1e2:
        print("WARNING: high basis correlation; K2/K3 may be sensitive to the fitting window.")
    print(f"Saved: {detail_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Fit a Bloch-like domain-wall profile.

Models
------
Mz(x)  = C - A tanh((x-x0)/delta)
Mxy(x) = C + A / cosh((x-x0)/delta)

The conventional Bloch-wall width is Delta_DW = pi * delta.
The script fits Mz and Mxy independently and reports both widths and their mean.

Input
-----
CSV or XLSX containing columns:
    x_minus_x0_nm, Mz, Mxy

This routine is intended only for Bloch-like walls. The broad low-temperature
2*theta0 / theta0 -> -theta0 metastable state should not be assigned a
conventional Bloch width with this model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported input format: {suffix}")


def mz_model(x, c, a, delta, x0):
    return c - a * np.tanh((x - x0) / delta)


def mxy_model(x, c, a, delta, x0):
    return c + a / np.cosh((x - x0) / delta)


def r_squared(y, y_fit):
    ss_res = np.sum((y - y_fit) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile", type=Path, help="Domain-wall profile CSV/XLSX")
    parser.add_argument("--x-column", default="x_minus_x0_nm")
    parser.add_argument("--mz-column", default="Mz")
    parser.add_argument("--mxy-column", default="Mxy")
    parser.add_argument("--output-prefix", default="bloch_fit")
    args = parser.parse_args()

    df = read_table(args.profile)
    required = [args.x_column, args.mz_column, args.mxy_column]
    missing = [name for name in required if name not in df.columns]
    if missing:
        raise KeyError(f"Missing required columns: {missing}")

    x = df[args.x_column].to_numpy(dtype=float)
    mz = df[args.mz_column].to_numpy(dtype=float)
    mxy = df[args.mxy_column].to_numpy(dtype=float)

    mask = np.isfinite(x) & np.isfinite(mz) & np.isfinite(mxy)
    x, mz, mxy = x[mask], mz[mask], mxy[mask]
    order = np.argsort(x)
    x, mz, mxy = x[order], mz[order], mxy[order]

    # Initial guesses. Offset terms are retained because finite-temperature
    # plane averages need not be perfectly symmetric around zero.
    p0_mz = [0.5 * (mz.max() + mz.min()), 0.5 * (mz.max() - mz.min()), 1.0, 0.0]
    popt_mz, pcov_mz = curve_fit(mz_model, x, mz, p0=p0_mz, maxfev=100000)

    delta_guess = max(abs(popt_mz[2]), 1e-6)
    p0_mxy = [mxy.min(), mxy.max() - mxy.min(), delta_guess, popt_mz[3]]
    popt_mxy, pcov_mxy = curve_fit(mxy_model, x, mxy, p0=p0_mxy, maxfev=100000)

    mz_fit = mz_model(x, *popt_mz)
    mxy_fit = mxy_model(x, *popt_mxy)

    delta_mz = abs(float(popt_mz[2]))
    delta_mxy = abs(float(popt_mxy[2]))
    width_mz = np.pi * delta_mz
    width_mxy = np.pi * delta_mxy
    width_mean = 0.5 * (width_mz + width_mxy)

    prefix = Path(args.output_prefix)
    detail_path = prefix.with_name(prefix.name + "_profile.csv")
    summary_path = prefix.with_name(prefix.name + "_summary.csv")

    pd.DataFrame({
        "x_minus_x0_nm": x,
        "Mz_data": mz,
        "Mz_fit": mz_fit,
        "Mxy_data": mxy,
        "Mxy_fit": mxy_fit,
    }).to_csv(detail_path, index=False)

    pd.DataFrame([{
        "input_file": str(args.profile),
        "Mz_C": popt_mz[0],
        "Mz_A": popt_mz[1],
        "Mz_delta_nm": delta_mz,
        "Mz_x0_nm": popt_mz[3],
        "Mz_R2": r_squared(mz, mz_fit),
        "Mxy_C": popt_mxy[0],
        "Mxy_A": popt_mxy[1],
        "Mxy_delta_nm": delta_mxy,
        "Mxy_x0_nm": popt_mxy[3],
        "Mxy_R2": r_squared(mxy, mxy_fit),
        "DW_width_from_Mz_nm": width_mz,
        "DW_width_from_Mxy_nm": width_mxy,
        "DW_width_mean_nm": width_mean,
    }]).to_csv(summary_path, index=False)

    print(f"Mz width   : {width_mz:.6f} nm")
    print(f"Mxy width  : {width_mxy:.6f} nm")
    print(f"Mean width : {width_mean:.6f} nm")
    print(f"Saved: {detail_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()

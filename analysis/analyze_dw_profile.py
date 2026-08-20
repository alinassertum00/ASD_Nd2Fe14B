#!/usr/bin/env python3
"""Post-process VAMPIRE domain-wall output into a 1-D averaged profile.

Expected VAMPIRE `dw-<step>.txt` layout used in this project:
    column 0 : x-bin index
    column 2 : Sx
    column 3 : Sy
    column 4 : Sz
    column 5 : statistical weight (if present)

For each x bin the script computes weighted averages Mx, My and Mz, then
Mxy=sqrt(Mx^2+My^2) and M=sqrt(Mx^2+My^2+Mz^2). The physical bin spacing is
read from `dimensions:system-size-x` in the VAMPIRE input and calculated as
Lx / number_of_x_bins. The slope-derived wall width is only a diagnostic; the
reported thesis width comes from the nonlinear Bloch fit in fit_bloch_wall.py.
"""

from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


def system_size_x_nm(input_file: Path):
    text = input_file.read_text()
    match = re.search(r"dimensions:system-size-x\s*=\s*([0-9.eE+-]+)", text)
    return float(match.group(1)) if match else None


def latest_dw_file(directory: Path):
    files = []
    for name in glob.glob(str(directory / "dw-*.txt")):
        match = re.search(r"dw-(\d+)\.txt$", name)
        if match:
            files.append((int(match.group(1)), Path(name)))
    if not files:
        raise FileNotFoundError(f"No dw-*.txt files found in {directory}")
    return max(files)[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dw-file", type=Path, default=None,
                   help="VAMPIRE dw-<step>.txt; default = latest in current directory")
    p.add_argument("--input-file", type=Path, default=Path("input"))
    p.add_argument("--left-fraction", type=float, default=0.15)
    p.add_argument("--right-fraction", type=float, default=0.15)
    p.add_argument("--output-dir", type=Path, default=Path("."))
    args = p.parse_args()

    dw_file = args.dw_file or latest_dw_file(Path("."))
    match = re.search(r"dw-(\d+)\.txt$", dw_file.name)
    step = int(match.group(1)) if match else -1

    data = np.loadtxt(dw_file)
    if data.ndim != 2 or data.shape[1] < 5:
        raise ValueError("Domain-wall output must contain at least five columns.")

    xbin = data[:, 0].astype(int)
    mx, my, mz = data[:, 2], data[:, 3], data[:, 4]
    weights = data[:, 5] if data.shape[1] > 5 else np.ones(len(data))
    xbins = np.unique(xbin)

    Mx, My, Mz = [], [], []
    for xb in xbins:
        mask = xbin == xb
        Mx.append(np.average(mx[mask], weights=weights[mask]))
        My.append(np.average(my[mask], weights=weights[mask]))
        Mz.append(np.average(mz[mask], weights=weights[mask]))

    Mx, My, Mz = map(np.asarray, (Mx, My, Mz))
    Mxy = np.sqrt(Mx**2 + My**2)
    M = np.sqrt(Mx**2 + My**2 + Mz**2)

    lx_nm = system_size_x_nm(args.input_file) if args.input_file.exists() else None
    if lx_nm is not None:
        dx_nm = lx_nm / len(xbins)
        x_nm = (xbins - xbins.min()) * dx_nm
    else:
        dx_nm = 1.0
        x_nm = xbins.astype(float)

    n = len(xbins)
    nl = max(1, int(args.left_fraction * n))
    nr = max(1, int(args.right_fraction * n))
    left_vec = np.array([Mx[:nl].mean(), My[:nl].mean(), Mz[:nl].mean()])
    right_vec = np.array([Mx[-nr:].mean(), My[-nr:].mean(), Mz[-nr:].mean()])

    left_unit = left_vec / np.linalg.norm(left_vec)
    right_unit = right_vec / np.linalg.norm(right_vec)
    domain_angle = math.degrees(math.acos(np.clip(np.dot(left_unit, right_unit), -1.0, 1.0)))

    # Standard/Bloch-like wall: Mz changes sign. Case-A-like wall: use Mx.
    if np.sign(left_vec[2]) != np.sign(right_vec[2]):
        wall_name, wall_var = "Mz", Mz
    else:
        wall_name, wall_var = "Mx", Mx

    left_value = wall_var[:nl].mean()
    right_value = wall_var[-nr:].mean()
    midpoint = 0.5 * (left_value + right_value)
    centre_index = int(np.argmin(np.abs(wall_var - midpoint)))
    x0_nm = float(x_nm[centre_index])

    gradient = np.gradient(wall_var, x_nm)
    max_slope = np.max(np.abs(gradient))
    slope_width = (abs(right_value - left_value) / max_slope if max_slope > 0 else np.nan)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    profile_path = args.output_dir / f"dw_profile_step_{step}.csv"
    summary_path = args.output_dir / f"dw_summary_step_{step}.csv"

    pd.DataFrame({
        "x_bin": xbins,
        "x_nm": x_nm,
        "x_minus_x0_nm": x_nm - x0_nm,
        "Mx": Mx,
        "My": My,
        "Mz": Mz,
        "Mxy": Mxy,
        "M": M,
    }).to_csv(profile_path, index=False)

    pd.DataFrame([{
        "file": str(dw_file),
        "step": step,
        "Lx_nm": lx_nm,
        "dx_nm": dx_nm,
        "x0_nm": x0_nm,
        "wall_variable": wall_name,
        "domain_wall_width_nm_slope_estimate": slope_width,
        "left_Mx": left_vec[0],
        "left_My": left_vec[1],
        "left_Mz": left_vec[2],
        "right_Mx": right_vec[0],
        "right_My": right_vec[1],
        "right_Mz": right_vec[2],
        "angle_between_domains_deg": domain_angle,
    }]).to_csv(summary_path, index=False)

    print(f"dx = {dx_nm:.6f} nm")
    print(f"x0 = {x0_nm:.6f} nm")
    print(f"wall variable = {wall_name}")
    print(f"domain angle = {domain_angle:.6f} deg")
    print(f"slope width diagnostic = {slope_width:.6f} nm")
    print(f"Saved: {profile_path}")
    print(f"Saved: {summary_path}")


if __name__ == "__main__":
    main()

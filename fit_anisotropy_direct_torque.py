#!/usr/bin/env python3
"""
Extract K1, K2, and K3 from constrained Monte Carlo torque data.

The fitted anisotropy energy is

    F(theta) = K1 sin^2(theta)
             + K2 sin^4(theta)
             + K3 sin^6(theta)

and the torque model is

    tau_y(theta) = -dF/dtheta

                   = -2 K1 sin(theta) cos(theta)
                     -4 K2 sin^3(theta) cos(theta)
                     -6 K3 sin^5(theta) cos(theta)

Required Excel columns
----------------------
Temperature
phi_deg
M
Ty_J

The physical rotation angle is phi_deg because theta_deg was fixed at zero
in the simulations.

Units
-----
Ty_J          : J
simulation V  : m^3
tau_y         : MJ/m^3
K1, K2, K3    : MJ/m^3
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# =====================================================================
# Frozen analysis settings
# =====================================================================

VOLUME_M3 = 1.1e-24

# Contiguous angular branches used in the final reproduction analysis.
#
# temperature: (lower angle, upper angle)
#
# These windows are explicit analysis inputs and should be reported
# transparently when presenting the fitted coefficients.
FIT_WINDOWS: Dict[int, Tuple[float, float]] = {
    25:  (3.0, 72.0),
    50:  (3.0, 66.0),
    75:  (33.0, 60.0),
    100: (15.0, 60.0),
    150: (3.0, 42.0),
    200: (30.0, 45.0),
    250: (18.0, 42.0),
    300: (30.0, 60.0),
    350: (33.0, 63.0),
    400: (30.0, 63.0),
    450: (27.0, 48.0),
    500: (36.0, 66.0),
}

REQUIRED_COLUMNS = {
    "Temperature",
    "phi_deg",
    "M",
    "Ty_J",
}


# =====================================================================
# Input handling
# =====================================================================

def temperature_from_filename(path: Path) -> int | None:
    """
    Extract a temperature from a filename such as:

        collected_150K_step3_10k.xlsx
    """
    match = re.search(r"(\d+)\s*K", path.name, flags=re.IGNORECASE)

    if match is None:
        return None

    return int(match.group(1))


def read_dataset(path: Path) -> pd.DataFrame:
    """Read and validate one collected CMC Excel file."""

    try:
        df = pd.read_excel(path)
    except Exception as exc:
        raise RuntimeError(
            f"Could not read Excel file {path}: {exc}"
        ) from exc

    missing = REQUIRED_COLUMNS.difference(df.columns)

    if missing:
        raise ValueError(
            f"{path.name} is missing required columns: "
            f"{sorted(missing)}"
        )

    df = df.copy()

    for column in REQUIRED_COLUMNS:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    invalid_rows = df[
        list(REQUIRED_COLUMNS)
    ].isna().any(axis=1)

    if invalid_rows.any():
        raise ValueError(
            f"{path.name} contains invalid numeric rows:\n"
            f"{df.loc[invalid_rows].to_string(index=False)}"
        )

    df = (
        df.sort_values("phi_deg")
        .drop_duplicates(
            subset="phi_deg",
            keep="last",
        )
        .reset_index(drop=True)
    )

    return df


# =====================================================================
# Fitting model
# =====================================================================

def torque_design_matrix(
    theta_rad: np.ndarray,
) -> np.ndarray:
    """
    Return the linear design matrix X for

        tau_y = X @ [K1, K2, K3].
    """

    sin_theta = np.sin(theta_rad)
    cos_theta = np.cos(theta_rad)

    basis_K1 = (
        -2.0
        * sin_theta
        * cos_theta
    )

    basis_K2 = (
        -4.0
        * sin_theta**3
        * cos_theta
    )

    basis_K3 = (
        -6.0
        * sin_theta**5
        * cos_theta
    )

    return np.column_stack(
        (
            basis_K1,
            basis_K2,
            basis_K3,
        )
    )


def fit_one_temperature(
    df: pd.DataFrame,
    temperature_K: int,
    lower_angle_deg: float,
    upper_angle_deg: float,
    volume_m3: float,
) -> tuple[dict, pd.DataFrame]:
    """
    Fit K1, K2, and K3 for one temperature.
    """

    selected = df[
        (df["phi_deg"] >= lower_angle_deg)
        & (df["phi_deg"] <= upper_angle_deg)
    ].copy()

    if len(selected) < 3:
        raise ValueError(
            f"{temperature_K} K has only "
            f"{len(selected)} points between "
            f"{lower_angle_deg} and "
            f"{upper_angle_deg} degrees."
        )

    theta_rad = np.deg2rad(
        selected["phi_deg"].to_numpy(
            dtype=float
        )
    )

    # Convert total torque in J into torque density in MJ/m^3.
    torque_density = (
        selected["Ty_J"].to_numpy(
            dtype=float
        )
        / volume_m3
        * 1.0e-6
    )

    X = torque_design_matrix(
        theta_rad
    )

    coefficients, _, rank, singular_values = (
        np.linalg.lstsq(
            X,
            torque_density,
            rcond=None,
        )
    )

    K1, K2, K3 = coefficients

    fitted_torque = (
        X @ coefficients
    )

    residuals = (
        torque_density
        - fitted_torque
    )

    number_of_points = len(
        torque_density
    )

    number_of_parameters = 3

    degrees_of_freedom = (
        number_of_points
        - number_of_parameters
    )

    sum_squared_error = float(
        np.sum(
            residuals**2
        )
    )

    rmse = float(
        np.sqrt(
            np.mean(
                residuals**2
            )
        )
    )

    mae = float(
        np.mean(
            np.abs(
                residuals
            )
        )
    )

    centered_torque = (
        torque_density
        - np.mean(
            torque_density
        )
    )

    total_sum_squares = float(
        np.sum(
            centered_torque**2
        )
    )

    if total_sum_squares > 0.0:
        r_squared = (
            1.0
            - sum_squared_error
            / total_sum_squares
        )
    else:
        r_squared = np.nan

    condition_number = float(
        np.linalg.cond(
            X
        )
    )

    # Classical least-squares covariance estimate.
    if (
        degrees_of_freedom > 0
        and rank == number_of_parameters
    ):
        residual_variance = (
            sum_squared_error
            / degrees_of_freedom
        )

        covariance = (
            residual_variance
            * np.linalg.inv(
                X.T @ X
            )
        )

        standard_errors = np.sqrt(
            np.diag(
                covariance
            )
        )

        confidence_intervals_95 = (
            1.96
            * standard_errors
        )

    else:
        standard_errors = np.full(
            3,
            np.nan,
        )

        confidence_intervals_95 = np.full(
            3,
            np.nan,
        )

    result = {
        "Temperature_K":
            temperature_K,

        "Lower_angle_deg":
            lower_angle_deg,

        "Upper_angle_deg":
            upper_angle_deg,

        "N_points":
            number_of_points,

        "K1_MJm3":
            float(K1),

        "K2_MJm3":
            float(K2),

        "K3_MJm3":
            float(K3),

        "K1_standard_error":
            float(
                standard_errors[0]
            ),

        "K2_standard_error":
            float(
                standard_errors[1]
            ),

        "K3_standard_error":
            float(
                standard_errors[2]
            ),

        "K1_CI95_half_width":
            float(
                confidence_intervals_95[0]
            ),

        "K2_CI95_half_width":
            float(
                confidence_intervals_95[1]
            ),

        "K3_CI95_half_width":
            float(
                confidence_intervals_95[2]
            ),

        "RMSE_MJm3":
            rmse,

        "MAE_MJm3":
            mae,

        "R_squared":
            float(
                r_squared
            ),

        "Condition_number":
            condition_number,

        "Minimum_M":
            float(
                selected["M"].min()
            ),

        "Maximum_M":
            float(
                selected["M"].max()
            ),

        "Mean_M":
            float(
                selected["M"].mean()
            ),
    }

    detailed_results = selected[
        [
            "phi_deg",
            "M",
            "Ty_J",
        ]
    ].copy()

    detailed_results.insert(
        0,
        "Temperature_K",
        temperature_K,
    )

    detailed_results[
        "theta_rad"
    ] = theta_rad

    detailed_results[
        "Ty_density_MJm3"
    ] = torque_density

    detailed_results[
        "Ty_fit_MJm3"
    ] = fitted_torque

    detailed_results[
        "residual_MJm3"
    ] = residuals

    return (
        result,
        detailed_results,
    )


# =====================================================================
# Plotting
# =====================================================================

def plot_torque_fit(
    detailed_results: pd.DataFrame,
    result: dict,
    output_directory: Path,
) -> None:
    """Save the torque data and fitted curve for one temperature."""

    temperature = int(
        result["Temperature_K"]
    )

    figure, axes = plt.subplots(
        figsize=(7.2, 4.8)
    )

    axes.scatter(
        detailed_results["phi_deg"],
        detailed_results[
            "Ty_density_MJm3"
        ],
        label="CMC torque",
        zorder=3,
    )

    axes.plot(
        detailed_results["phi_deg"],
        detailed_results[
            "Ty_fit_MJm3"
        ],
        label="Direct anisotropy fit",
        linewidth=1.8,
    )

    axes.set_xlabel(
        r"Rotation angle $\phi$ (degree)"
    )

    axes.set_ylabel(
        r"Torque density "
        r"$\tau_y$ "
        r"(MJ m$^{-3}$)"
    )

    axes.set_title(
        f"{temperature} K: "
        f"{result['Lower_angle_deg']:g}°"
        f"–"
        f"{result['Upper_angle_deg']:g}°"
    )

    axes.grid(
        True,
        alpha=0.25,
    )

    axes.legend()

    figure.tight_layout()

    figure.savefig(
        output_directory
        / f"torque_fit_{temperature}K.png",
        dpi=300,
    )

    plt.close(
        figure
    )


# =====================================================================
# File discovery
# =====================================================================

def discover_input_files(
    input_directory: Path,
    pattern: str,
) -> list[Path]:
    """Find and sort all input Excel files."""

    files = list(
        input_directory.glob(
            pattern
        )
    )

    if not files:
        raise FileNotFoundError(
            f"No files matching "
            f"'{pattern}' were found in "
            f"{input_directory}."
        )

    def sorting_key(
        path: Path,
    ) -> tuple[int, str]:

        temperature = (
            temperature_from_filename(
                path
            )
        )

        if temperature is None:
            temperature = 10**9

        return (
            temperature,
            path.name,
        )

    return sorted(
        files,
        key=sorting_key,
    )


# =====================================================================
# Main program
# =====================================================================

def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Fit K1, K2, and K3 "
            "from CMC Ty torque data."
        )
    )

    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("."),
        help=(
            "Directory containing "
            "the collected Excel files."
        ),
    )

    parser.add_argument(
        "--pattern",
        default="collected_*K*.xlsx",
        help=(
            "Glob pattern used to locate "
            "the input Excel files."
        ),
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(
            "anisotropy_fit_results"
        ),
        help=(
            "Directory in which the "
            "results will be saved."
        ),
    )

    parser.add_argument(
        "--volume",
        type=float,
        default=VOLUME_M3,
        help=(
            "Simulation volume in m^3."
        ),
    )

    parser.add_argument(
        "--no-plots",
        action="store_true",
        help=(
            "Disable PNG plot generation."
        ),
    )

    arguments = parser.parse_args()

    arguments.output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    input_files = (
        discover_input_files(
            arguments.input_dir,
            arguments.pattern,
        )
    )

    all_summary_results = []
    all_detailed_results = []

    print(
        "Direct torque fitting "
        "of K1, K2, and K3"
    )

    print(
        f"Simulation volume = "
        f"{arguments.volume:.6e} m^3"
    )

    print()

    for input_file in input_files:

        dataframe = read_dataset(
            input_file
        )

        temperature_from_data = int(
            round(
                float(
                    dataframe[
                        "Temperature"
                    ].median()
                )
            )
        )

        temperature_from_name = (
            temperature_from_filename(
                input_file
            )
        )

        if (
            temperature_from_name
            is not None
            and temperature_from_name
            != temperature_from_data
        ):
            raise ValueError(
                f"{input_file.name}: "
                f"filename indicates "
                f"{temperature_from_name} K, "
                f"but the spreadsheet "
                f"contains "
                f"{temperature_from_data} K."
            )

        temperature = (
            temperature_from_data
        )

        if temperature not in FIT_WINDOWS:
            print(
                f"Skipping "
                f"{input_file.name}: "
                f"no frozen fitting window "
                f"for {temperature} K."
            )

            continue

        lower_angle, upper_angle = (
            FIT_WINDOWS[
                temperature
            ]
        )

        result, detailed_result = (
            fit_one_temperature(
                df=dataframe,
                temperature_K=temperature,
                lower_angle_deg=lower_angle,
                upper_angle_deg=upper_angle,
                volume_m3=arguments.volume,
            )
        )

        result[
            "Input_file"
        ] = input_file.name

        all_summary_results.append(
            result
        )

        all_detailed_results.append(
            detailed_result
        )

        if not arguments.no_plots:
            plot_torque_fit(
                detailed_result,
                result,
                arguments.output_dir,
            )

        print(
            f"{temperature:>3d} K | "
            f"{lower_angle:>4.0f}°–"
            f"{upper_angle:<4.0f}° | "
            f"K1 = "
            f"{result['K1_MJm3']:>8.3f} | "
            f"K2 = "
            f"{result['K2_MJm3']:>8.3f} | "
            f"K3 = "
            f"{result['K3_MJm3']:>8.3f} | "
            f"R² = "
            f"{result['R_squared']:.6f}"
        )

    if not all_summary_results:
        raise RuntimeError(
            "No datasets were fitted."
        )

    summary_dataframe = pd.DataFrame(
        all_summary_results
    ).sort_values(
        "Temperature_K"
    )

    detailed_dataframe = pd.concat(
        all_detailed_results,
        ignore_index=True,
    )

    summary_csv = (
        arguments.output_dir
        / "anisotropy_K_summary.csv"
    )

    details_csv = (
        arguments.output_dir
        / "anisotropy_fit_points.csv"
    )

    excel_output = (
        arguments.output_dir
        / "anisotropy_fit_results.xlsx"
    )

    summary_dataframe.to_csv(
        summary_csv,
        index=False,
    )

    detailed_dataframe.to_csv(
        details_csv,
        index=False,
    )

    with pd.ExcelWriter(
        excel_output,
        engine="openpyxl",
    ) as writer:

        summary_dataframe.to_excel(
            writer,
            sheet_name="K_summary",
            index=False,
        )

        detailed_dataframe.to_excel(
            writer,
            sheet_name="fit_points",
            index=False,
        )

        for (
            temperature,
            temperature_data,
        ) in detailed_dataframe.groupby(
            "Temperature_K"
        ):

            sheet_name = (
                f"{int(temperature)}K"
            )

            temperature_data.to_excel(
                writer,
                sheet_name=sheet_name,
                index=False,
            )

    print()

    print("Saved:")

    print(
        f"  {summary_csv}"
    )

    print(
        f"  {details_csv}"
    )

    print(
        f"  {excel_output}"
    )

    return 0


if __name__ == "__main__":

    try:
        raise SystemExit(
            main()
        )

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        raise SystemExit(1)

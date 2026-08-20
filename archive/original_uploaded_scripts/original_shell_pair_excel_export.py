import numpy as np
import pandas as pd
from itertools import product
from pathlib import Path
from math import ceil

UCF_FILE = "Nd2Fe14B_scaled_175.ucf"
CUTOFF_RADIUS_A = 3.376376945
SHELL_TOL_A = 1e-4
OUTPUT_EXCEL = "Nd2Fe14B_shell_pairs.xlsx"

SITE_MAP = {
    0: "Nd(4g)",
    1: "Nd(4f)",
    2: "Fe(4c)",
    3: "Fe(4e)",
    4: "Fe(8j2)",
    5: "Fe(8j1)",
    6: "Fe(16k1)",
    7: "Fe(16k2)",
    8: "B"
}

FE_MATS = {2, 3, 4, 5, 6, 7}
ND_MATS = {0, 1}

def pair_type(mi, mj):
    if mi in FE_MATS and mj in FE_MATS:
        return "Fe-Fe"
    if (mi in FE_MATS and mj in ND_MATS) or (mi in ND_MATS and mj in FE_MATS):
        return "Fe-Nd"
    if mi in ND_MATS and mj in ND_MATS:
        return "Nd-Nd"
    if mi == 8 or mj == 8:
        return "B-related"
    return "Other"

def keep_pair(mi, mj):
    return pair_type(mi, mj) in {"Fe-Fe", "Fe-Nd"}

def read_ucf(path):
    lines = Path(path).read_text().splitlines()

    cell = None
    for i, line in enumerate(lines):
        if "Unit cell size" in line:
            cell = np.array([float(x) for x in lines[i+1].split()[:3]])
            break

    if cell is None:
        raise RuntimeError("Could not find unit-cell size in UCF.")

    atom_start = None
    num_atoms = None

    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
            possible_n = int(parts[0])
            if possible_n > 10:
                atom_start = i + 1
                num_atoms = possible_n
                break

    if atom_start is None:
        raise RuntimeError("Could not find atom section in UCF.")

    atoms = []
    for line in lines[atom_start:atom_start + num_atoms]:
        p = line.split()
        atoms.append({
            "id": int(p[0]),
            "frac": np.array([float(p[1]), float(p[2]), float(p[3])]),
            "mat": int(p[4]),
            "cat": int(p[5]),
            "hcat": int(p[6]),
            "site": SITE_MAP.get(int(p[4]), f"mat{p[4]}")
        })

    return cell, atoms

def assign_shells(distances, tol):
    shell_centres = []
    shell_ids = []

    for d in distances:
        assigned = False
        for k, c in enumerate(shell_centres):
            if abs(d - c) <= tol:
                shell_ids.append(k + 1)
                assigned = True
                break
        if not assigned:
            shell_centres.append(d)
            shell_ids.append(len(shell_centres))

    return shell_ids

def main():
    if not Path(UCF_FILE).exists():
        raise FileNotFoundError(
            f"Could not find {UCF_FILE}. Put it in this folder or edit UCF_FILE."
        )

    cell, atoms = read_ucf(UCF_FILE)

    max_image = int(ceil(CUTOFF_RADIUS_A / min(cell))) + 1
    image_range = range(-max_image, max_image + 1)

    unique = {}

    for ai in atoms:
        for aj in atoms:
            if not keep_pair(ai["mat"], aj["mat"]):
                continue

            for dx, dy, dz in product(image_range, repeat=3):
                if ai["id"] == aj["id"] and (dx, dy, dz) == (0, 0, 0):
                    continue

                shift = np.array([dx, dy, dz])
                dr_frac = aj["frac"] + shift - ai["frac"]
                dr_cart = dr_frac * cell
                dist = np.linalg.norm(dr_cart)

                if dist <= CUTOFF_RADIUS_A + 1e-10:
                    key1 = (ai["id"], aj["id"], dx, dy, dz)
                    key2 = (aj["id"], ai["id"], -dx, -dy, -dz)
                    canonical = min(key1, key2)

                    if canonical not in unique:
                        unique[canonical] = {
                            "i": ai["id"],
                            "j": aj["id"],
                            "dx": dx,
                            "dy": dy,
                            "dz": dz,
                            "distance_A": dist,
                            "mat_i": ai["mat"],
                            "mat_j": aj["mat"],
                            "site_i": ai["site"],
                            "site_j": aj["site"],
                            "pair_type": pair_type(ai["mat"], aj["mat"])
                        }

    df = pd.DataFrame(unique.values())
    df = df.sort_values("distance_A").reset_index(drop=True)

    df["shell"] = assign_shells(df["distance_A"].values, SHELL_TOL_A)
    df["site_pair"] = df.apply(
        lambda r: "--".join(sorted([r["site_i"], r["site_j"]])),
        axis=1
    )

    directed_rows = []
    for _, r in df.iterrows():
        directed_rows.append(r.to_dict())

        rr = r.to_dict()
        rr["i"], rr["j"] = r["j"], r["i"]
        rr["mat_i"], rr["mat_j"] = r["mat_j"], r["mat_i"]
        rr["site_i"], rr["site_j"] = r["site_j"], r["site_i"]
        rr["dx"], rr["dy"], rr["dz"] = -r["dx"], -r["dy"], -r["dz"]
        directed_rows.append(rr)

    df_directed = pd.DataFrame(directed_rows)

    shell_summary = (
        df.groupby(["shell", "pair_type"])
        .size()
        .reset_index(name="unique_pairs")
    )
    shell_summary["directed_lines"] = 2 * shell_summary["unique_pairs"]

    site_summary = (
        df.groupby(["shell", "pair_type", "site_pair"])
        .size()
        .reset_index(name="unique_pairs")
    )
    site_summary["directed_lines"] = 2 * site_summary["unique_pairs"]

    overall = pd.DataFrame([{
        "ucf_file": UCF_FILE,
        "cutoff_radius_A": CUTOFF_RADIUS_A,
        "shell_tolerance_A": SHELL_TOL_A,
        "unique_pairs": len(df),
        "symmetrized_directed_lines": len(df_directed),
        "max_distance_A": df["distance_A"].max(),
        "number_of_shells_found": df["shell"].max()
    }])

    with pd.ExcelWriter(OUTPUT_EXCEL, engine="openpyxl") as writer:
        overall.to_excel(writer, sheet_name="overall", index=False)
        shell_summary.to_excel(writer, sheet_name="shell_summary", index=False)
        site_summary.to_excel(writer, sheet_name="site_pair_summary", index=False)
        df.to_excel(writer, sheet_name="unique_pairs", index=False)
        df_directed.to_excel(writer, sheet_name="directed_ucf_lines", index=False)

    print("DONE")
    print(f"Cutoff radius: {CUTOFF_RADIUS_A:.9f} A")
    print(f"Unique pairs: {len(df)}")
    print(f"Directed UCF lines: {len(df_directed)}")
    print(f"Excel written to: {OUTPUT_EXCEL}")

if __name__ == "__main__":
    main()

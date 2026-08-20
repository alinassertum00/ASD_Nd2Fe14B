#!/usr/bin/env python3
"""Reconstruct selected Nd2Fe14B neighbour pairs from a VAMPIRE UCF file.

The script is portable: all file paths and numerical settings are CLI options.
It keeps Fe-Fe and Fe-Nd pairs by default and exports an Excel workbook with
shell, site-pair and directed-line summaries.
"""
from __future__ import annotations
import argparse
from itertools import product
from math import ceil
from pathlib import Path
import numpy as np
import pandas as pd

SITE_MAP = {0:"Nd(4g)",1:"Nd(4f)",2:"Fe(4c)",3:"Fe(4e)",4:"Fe(8j2)",5:"Fe(8j1)",6:"Fe(16k1)",7:"Fe(16k2)",8:"B"}
FE_MATS = {2,3,4,5,6,7}
ND_MATS = {0,1}

def pair_type(mi, mj):
    if mi in FE_MATS and mj in FE_MATS: return "Fe-Fe"
    if (mi in FE_MATS and mj in ND_MATS) or (mi in ND_MATS and mj in FE_MATS): return "Fe-Nd"
    if mi in ND_MATS and mj in ND_MATS: return "Nd-Nd"
    if mi == 8 or mj == 8: return "B-related"
    return "Other"

def read_ucf(path: Path):
    lines = path.read_text().splitlines()
    cell = None
    for i, line in enumerate(lines):
        if "Unit cell size" in line:
            cell = np.array([float(x) for x in lines[i+1].split()[:3]])
            break
    if cell is None:
        raise RuntimeError("Could not find unit-cell size in UCF.")
    atom_start = num_atoms = None
    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit() and int(parts[0]) > 10:
            atom_start, num_atoms = i + 1, int(parts[0]); break
    if atom_start is None:
        raise RuntimeError("Could not find atom section in UCF.")
    atoms = []
    for line in lines[atom_start:atom_start+num_atoms]:
        p = line.split(); mat = int(p[4])
        atoms.append({"id":int(p[0]),"frac":np.array([float(p[1]),float(p[2]),float(p[3])]),"mat":mat,"site":SITE_MAP.get(mat,f"mat{mat}")})
    return cell, atoms

def assign_shells(distances, tol):
    centres, ids = [], []
    for d in distances:
        for k, c in enumerate(centres):
            if abs(d-c) <= tol:
                ids.append(k+1); break
        else:
            centres.append(d); ids.append(len(centres))
    return ids

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("ucf", type=Path, help="Input VAMPIRE UCF file")
    p.add_argument("--cutoff-A", type=float, default=3.376376945)
    p.add_argument("--shell-tol-A", type=float, default=1e-4)
    p.add_argument("--output", type=Path, default=Path("Nd2Fe14B_shell_pairs.xlsx"))
    args = p.parse_args()
    ucf = args.ucf.expanduser().resolve()
    cell, atoms = read_ucf(ucf)
    max_image = int(ceil(args.cutoff_A/min(cell))) + 1
    image_range = range(-max_image, max_image+1)
    unique = {}
    for ai in atoms:
        for aj in atoms:
            if pair_type(ai["mat"], aj["mat"]) not in {"Fe-Fe","Fe-Nd"}: continue
            for dx,dy,dz in product(image_range, repeat=3):
                if ai["id"] == aj["id"] and (dx,dy,dz)==(0,0,0): continue
                dr = (aj["frac"] + np.array([dx,dy,dz]) - ai["frac"]) * cell
                dist = np.linalg.norm(dr)
                if dist <= args.cutoff_A + 1e-10:
                    k1=(ai["id"],aj["id"],dx,dy,dz); k2=(aj["id"],ai["id"],-dx,-dy,-dz); key=min(k1,k2)
                    unique.setdefault(key,{"i":ai["id"],"j":aj["id"],"dx":dx,"dy":dy,"dz":dz,"distance_A":dist,"mat_i":ai["mat"],"mat_j":aj["mat"],"site_i":ai["site"],"site_j":aj["site"],"pair_type":pair_type(ai["mat"],aj["mat"])})
    df = pd.DataFrame(unique.values()).sort_values("distance_A").reset_index(drop=True)
    df["shell"] = assign_shells(df["distance_A"].values, args.shell_tol_A)
    df["site_pair"] = df.apply(lambda r:"--".join(sorted([r["site_i"],r["site_j"]])),axis=1)
    directed=[]
    for _,r in df.iterrows():
        directed.append(r.to_dict()); rr=r.to_dict(); rr["i"],rr["j"]=r["j"],r["i"]; rr["mat_i"],rr["mat_j"]=r["mat_j"],r["mat_i"]; rr["site_i"],rr["site_j"]=r["site_j"],r["site_i"]; rr["dx"],rr["dy"],rr["dz"]=-r["dx"],-r["dy"],-r["dz"]; directed.append(rr)
    df_directed=pd.DataFrame(directed)
    shell_summary=df.groupby(["shell","pair_type"]).size().reset_index(name="unique_pairs"); shell_summary["directed_lines"]=2*shell_summary["unique_pairs"]
    site_summary=df.groupby(["shell","pair_type","site_pair"]).size().reset_index(name="unique_pairs"); site_summary["directed_lines"]=2*site_summary["unique_pairs"]
    overall=pd.DataFrame([{"ucf_file":str(ucf),"cutoff_radius_A":args.cutoff_A,"shell_tolerance_A":args.shell_tol_A,"unique_pairs":len(df),"symmetrized_directed_lines":len(df_directed),"max_distance_A":df["distance_A"].max(),"number_of_shells_found":df["shell"].max()}])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(args.output, engine="openpyxl") as w:
        overall.to_excel(w,"overall",index=False); shell_summary.to_excel(w,"shell_summary",index=False); site_summary.to_excel(w,"site_pair_summary",index=False); df.to_excel(w,"unique_pairs",index=False); df_directed.to_excel(w,"directed_ucf_lines",index=False)
    print(f"Unique pairs: {len(df)}"); print(f"Directed UCF lines: {len(df_directed)}"); print(f"Saved: {args.output}")
if __name__ == "__main__": main()

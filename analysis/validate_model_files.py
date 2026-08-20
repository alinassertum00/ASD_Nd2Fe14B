#!/usr/bin/env python3

import argparse
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate basic consistency of VAMPIRE MAT and UCF files."
    )
    parser.add_argument("mat", type=Path, help="VAMPIRE material file")
    parser.add_argument("ucf", type=Path, help="VAMPIRE unit-cell file")
    return parser.parse_args()


def main():
    args = parse_args()

    mat = args.mat.expanduser()
    ucf = args.ucf.expanduser()

    if not mat.is_file():
        raise SystemExit(f"MAT file not found: {mat}")

    if not ucf.is_file():
        raise SystemExit(f"UCF file not found: {ucf}")

    mt = mat.read_text(encoding="utf-8", errors="replace")
    ut = ucf.read_text(encoding="utf-8", errors="replace")

    mat_material_match = re.search(
        r"material:num-materials\s*=\s*(\d+)",
        mt,
        flags=re.IGNORECASE,
    )
    mat_materials = (
        int(mat_material_match.group(1))
        if mat_material_match
        else None
    )

    exchange_match = re.search(
        r"(\d+)\s+isotropic",
        ut,
        flags=re.IGNORECASE,
    )
    exchange_count = (
        int(exchange_match.group(1))
        if exchange_match
        else None
    )

    has_rescaling = bool(
        re.search(
            r"temperature-rescaling",
            mt,
            flags=re.IGNORECASE,
        )
    )

    has_dw_initialization = bool(
        re.search(
            r"domain-wall-second-magnetisation-vector|initial-spin-direction",
            mt,
            flags=re.IGNORECASE,
        )
    )

    print(f"MAT materials declared: {mat_materials}")
    print(
        "UCF directed isotropic interactions declared: "
        f"{exchange_count}"
    )
    print(f"Contains temperature rescaling: {has_rescaling}")
    print(
        "Contains domain-wall initialization in MAT: "
        f"{has_dw_initialization}"
    )


if __name__ == "__main__":
    main()

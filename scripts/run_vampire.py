#!/usr/bin/env python3
"""Run one VAMPIRE simulation directory in a machine-independent way.

The VAMPIRE executable is resolved in this order:
1. --vampire command-line option
2. VAMPIRE_BIN environment variable
3. PATH lookup for vampire-serial

Example:
    python scripts/run_vampire.py simulations/domain_wall --vampire /opt/vampire/vampire-serial
"""
from __future__ import annotations
import argparse
import os
import shutil
import subprocess
from pathlib import Path


def resolve_executable(value: str | None) -> str:
    candidate = value or os.environ.get("VAMPIRE_BIN") or "vampire-serial"
    expanded = str(Path(candidate).expanduser())
    if Path(expanded).is_file():
        return expanded
    found = shutil.which(candidate)
    if found:
        return found
    raise FileNotFoundError(
        f"Cannot find VAMPIRE executable '{candidate}'. Pass --vampire, set VAMPIRE_BIN, "
        "or add vampire-serial to PATH."
    )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case_dir", type=Path, help="Directory containing a VAMPIRE input file")
    p.add_argument("--vampire", default=None, help="Path/name of VAMPIRE executable")
    p.add_argument("--log", default="run.log", help="Log filename inside case directory")
    args = p.parse_args()

    case_dir = args.case_dir.expanduser().resolve()
    if not (case_dir / "input").exists():
        raise FileNotFoundError(f"No 'input' file found in {case_dir}")

    exe = resolve_executable(args.vampire)
    log_path = case_dir / args.log
    print(f"VAMPIRE : {exe}")
    print(f"Case    : {case_dir}")
    print(f"Log     : {log_path}")
    with log_path.open("w") as log:
        subprocess.run([exe], cwd=case_dir, stdout=log, stderr=subprocess.STDOUT, check=True)
    print("Simulation completed successfully.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run all angle folders in a generated CMC sweep with bounded parallelism.

Portable replacement for VM-specific shell loops. The VAMPIRE executable is
resolved from --vampire, VAMPIRE_BIN, or PATH (vampire-serial).
"""
from __future__ import annotations
import argparse
import concurrent.futures
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


def run_case(exe: str, case: Path) -> tuple[str, int]:
    with (case / "run.log").open("w") as log:
        proc = subprocess.run([exe], cwd=case, stdout=log, stderr=subprocess.STDOUT)
    return case.name, proc.returncode


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("sweep_dir", type=Path)
    p.add_argument("--vampire", default=None)
    p.add_argument("--jobs", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    args = p.parse_args()

    sweep = args.sweep_dir.expanduser().resolve()
    cases = sorted(p for p in sweep.glob("angle_*deg") if (p / "input").exists())
    if not cases:
        raise FileNotFoundError(f"No angle_*deg/input cases found in {sweep}")
    exe = resolve_executable(args.vampire)
    print(f"Running {len(cases)} CMC cases with up to {args.jobs} concurrent jobs")
    print(f"VAMPIRE: {exe}")

    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_case, exe, c) for c in cases]
        for future in concurrent.futures.as_completed(futures):
            name, code = future.result()
            print(f"{name}: {'OK' if code == 0 else f'FAILED ({code})'}")
            if code != 0:
                failures.append(name)
    if failures:
        raise SystemExit(f"Failed cases: {', '.join(failures)}")
    print("All CMC angle jobs completed successfully.")


if __name__ == "__main__":
    main()

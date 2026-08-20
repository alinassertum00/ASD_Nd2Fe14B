#!/usr/bin/env python3
"""Generate VAMPIRE constrained-Monte-Carlo angle folders for one temperature.

The thesis workflow rotates the constrained magnetisation through phi while
keeping theta fixed at zero in VAMPIRE's coordinate convention. Each angle is
an independent simulation and can therefore be run in parallel.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--temperature", type=float, required=True)
    p.add_argument("--start-angle", type=float, default=0.0)
    p.add_argument("--stop-angle", type=float, default=90.0)
    p.add_argument("--step-angle", type=float, default=3.0)
    p.add_argument("--equilibration-steps", type=int, default=10000)
    p.add_argument("--averaging-steps", type=int, default=10000)
    p.add_argument("--system-size-nm", type=float, default=10.0)
    p.add_argument("--material", type=Path, default=PROJECT_ROOT / "model" / "Nd2Fe14B.mat")
    p.add_argument("--ucf", type=Path, default=PROJECT_ROOT / "model" / "Nd2Fe14B_scaled_175.ucf")
    p.add_argument("--output", type=Path, default=None)
    args = p.parse_args()

    output = args.output or Path(f"CMC_{args.temperature:g}K")
    output.mkdir(parents=True, exist_ok=True)
    shutil.copy2(args.material, output / "Nd2Fe14B.mat")
    shutil.copy2(args.ucf, output / "Nd2Fe14B_scaled_175.ucf")

    n = int(round((args.stop_angle - args.start_angle) / args.step_angle))
    angles = [args.start_angle + i * args.step_angle for i in range(n + 1)]

    for angle in angles:
        folder = output / f"angle_{angle:07.3f}deg"
        folder.mkdir(exist_ok=True)
        text = f"""material:unit-cell-file = ../Nd2Fe14B_scaled_175.ucf
material:file = ../Nd2Fe14B.mat

create:periodic-boundaries-x
create:periodic-boundaries-y
create:periodic-boundaries-z

dimensions:system-size-x = {args.system_size_nm:g} !nm
dimensions:system-size-y = {args.system_size_nm:g} !nm
dimensions:system-size-z = {args.system_size_nm:g} !nm

sim:minimum-temperature = {args.temperature:g}
sim:maximum-temperature = {args.temperature:g}
sim:temperature-increment = 1

sim:equilibration-time-steps = {args.equilibration_steps}
sim:loop-time-steps = {args.averaging_steps}

sim:constraint-angle-theta-minimum = 0
sim:constraint-angle-theta-maximum = 0
sim:constraint-angle-theta-increment = 1

sim:constraint-angle-phi-minimum = {angle:.6f}
sim:constraint-angle-phi-maximum = {angle:.6f}
sim:constraint-angle-phi-increment = 1

sim:integrator = constrained-monte-carlo
sim:program = cmc-anisotropy

output:temperature
output:constraint-theta
output:constraint-phi
output:mean-magnetisation-length
output:mean-total-torque
output:mean-total-energy
"""
        (folder / "input").write_text(text)

    print(f"Created {len(angles)} angle cases in {output}")


if __name__ == "__main__":
    main()

"""Parallel fixed-rule finite-difference columns for the canonical bulk AR4 charging gate."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import os
import subprocess
import sys


parser = argparse.ArgumentParser()
parser.add_argument("--checkpoint", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--step", type=float, required=True)
parser.add_argument("--workers", type=int, default=5)
parser.add_argument("--adjoint-max", type=int, default=20)
parser.add_argument("--forward-max", type=int, default=20)
parser.add_argument("--element-absolute", type=float, default=0.01)
parser.add_argument("--element-relative", type=float, default=0.15)
parser.add_argument(
    "--driver", default=str(Path(__file__).with_name("charging_nodal_campaign.py")))
args = parser.parse_args()
if (args.step == 0.0 or args.workers <= 0 or args.adjoint_max < 0
        or args.forward_max < 0 or args.element_absolute < 0.0
        or args.element_relative < 0.0):
    raise ValueError("campaign step, workers, quadrature levels, and tolerances are invalid")

output_directory = Path(args.output)
output_directory.mkdir(parents=True, exist_ok=True)
direction = "plus" if args.step > 0.0 else "minus"

base = [
    sys.executable, args.driver,
    "--geometry", "bulk", "--trench-width", "4", "--trench-depth", "16",
    "--side-thickness", "5", "--nodal", "--poisson", "--update", "picard",
    "--iterations", "1", "--beta", "0.025", "--override-restart-beta", "0.025",
    "--fixed-dt", "0.01", "--face-offset", "0.000001", "--grazing",
    "--freeze-method", "--freeze-levels",
    "--adjoint-max", str(args.adjoint_max), "--forward-max", str(args.forward_max),
    "--element-absolute", str(args.element_absolute),
    "--element-relative", str(args.element_relative),
    "--initial", args.checkpoint,
]
environment = dict(os.environ)
environment.setdefault("PETCH_DEVICE", "cuda")


def run(dof):
    output = output_directory / f"dof_{dof:02d}_{direction}.npz"
    log = output_directory / f"dof_{dof:02d}_{direction}.log"
    command = base + [
        "--perturb-charge-dof", str(dof),
        "--perturb-charge-coordinate-volts", str(args.step),
        "--output", str(output),
    ]
    completed = subprocess.run(
        command, env=environment, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, check=False)
    log.write_text(completed.stdout)
    if completed.returncode != 0:
        raise RuntimeError(f"dof {dof} failed; see {log}")
    return dof


with ThreadPoolExecutor(max_workers=args.workers) as pool:
    futures = {pool.submit(run, dof): dof for dof in range(47)}
    for future in as_completed(futures):
        dof = future.result()
        print(f"completed {dof}", flush=True)

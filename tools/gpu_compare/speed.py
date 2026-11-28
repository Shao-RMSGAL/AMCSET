#!/usr/bin/env python3
# Copyright 2025, Texas A&M University
#
# This file is part of AMCSET.
#
# AMCSET is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# AMCSET is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# AMCSET. If not, see <https://www.gnu.org/licenses/>.
"""Wall time of AMCSET on the CPU and on GPU backends (one run each).

Runs the cases of compare_backends.py at the sizes below with each backend,
and records the wall time of `amcset run` and the run's "compute" report
(summary.json: setup, kernel, transfer and host times, launches, steps) in
<out>/speed.json, which make_figures.py draws. Indications, not a benchmark:
time on an idle machine for numbers to quote (docs/GPU.md "Speed").

  python tools/gpu_compare/speed.py --amcset build/Release/amcset.exe --out build/gpu_compare/speed
"""

import argparse
import json
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_backends import CASES  # noqa: E402

ALL = ["cpu", "opencl:double", "opencl:mixed", "hip:double", "hip:mixed"]
MIXED = ["cpu", "opencl:mixed", "hip:mixed"]
# (case, ions, backends): double precision only where it finishes in about a minute.
RUNS = [
    ("fe50keV_kp", 1_000_000, ALL),
    ("fe50keV_full", 200_000, ALL),
    ("h100keV_Si_full", 500_000, ALL),
    ("ar30keV_SiO2_Si_full", 200_000, MIXED),
    ("au200keV_Si_full", 16_000, MIXED),
    ("fe50keV_full", 20_000, MIXED),
    ("fe1MeV_full", 2_000, MIXED),
]


def live_series(a, scenarios):
    """--live: the first N ions of one run for growing N (make_figures.py's live animation).

    Blocks of 16 ions are merged in block order and every ion has its own random
    stream, so a run of N ions with a seed is exactly what a longer run with that
    seed shows after its first N ions: these are the live plot's snapshots."""
    sizes = [int(1600 * 2 ** (k / 2)) // 16 * 16 for k in range(15)]
    name, precision = a.live.split(":") if a.live != "cpu" else ("cpu", None)
    for n in sizes:
        d = os.path.join(a.out, f"{n:07d}")
        os.makedirs(d, exist_ok=True)
        sc = dict(scenarios["fe50keV_full"])
        sc["stopping"] = {"model": "srim"}
        path = os.path.join(d, "scenario.json")
        with open(path, "w") as f:
            json.dump(sc, f)
        args = [a.amcset, "run", path, "-q", "-n", str(n), "-s", str(a.seed), "-o", d, "--formats", "trim",
                "--set", "output.exit_lists=false", "--set", "output.range_3d=false"]
        if precision:
            args += ["--compute", name, "--precision", precision]
        subprocess.run(args, check=True, capture_output=True)
        print(f"live {n:7d} ions", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amcset", required=True)
    ap.add_argument("--out", default="build/gpu_compare/speed")
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--live", metavar="BACKEND",
                    help="instead: runs of the first N ions of 50 keV Fe (full cascades) on BACKEND (e.g. hip:mixed)")
    a = ap.parse_args()
    scenarios = {name: sc for name, sc, _ in CASES}
    if a.live:
        live_series(a, scenarios)
        return 0
    results = []
    for case, ions, backends in RUNS:
        for b in backends:
            d = os.path.join(a.out, f"{case}_{ions}", b.replace(":", "_"))
            os.makedirs(d, exist_ok=True)
            sc = dict(scenarios[case])
            sc["stopping"] = {"model": "srim"}
            path = os.path.join(d, "scenario.json")
            with open(path, "w") as f:
                json.dump(sc, f)
            args = [a.amcset, "run", path, "-q", "-n", str(ions), "-s", str(a.seed), "-o", d, "--formats", "trim",
                    "--set", "output.exit_lists=false", "--set", "output.range_3d=false"]
            if b != "cpu":
                name, precision = b.split(":")
                args += ["--compute", name, "--precision", precision]
            t0 = time.time()
            p = subprocess.run(args, capture_output=True, text=True)
            wall = time.time() - t0
            if p.returncode != 0:
                print(f"{case} {ions} {b}: failed\n{p.stderr}", file=sys.stderr)
                continue
            with open(os.path.join(d, "summary.json")) as f:
                s = json.load(f)
            row = {"case": case, "ions": ions, "backend": b, "wall_seconds": wall,
                   "elapsed_seconds": s.get("elapsed_seconds"), "compute": s.get("compute")}
            results.append(row)
            print(f"{case:22s} {ions:8d} {b:14s} {wall:7.2f} s", flush=True)
            with open(os.path.join(a.out, "speed.json"), "w") as f:
                json.dump(results, f, indent=1)
    print(f"Written: {os.path.join(a.out, 'speed.json')}")


if __name__ == "__main__":
    sys.exit(main())

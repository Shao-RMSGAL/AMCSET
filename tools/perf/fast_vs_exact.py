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

"""Statistical check of the fast approximations (physics.fast_math).

Runs every reference case of bench.py with several seeds, once exactly and
once with --fast, and compares every number of summary.json:

    python tools/perf/fast_vs_exact.py --exe build/Release/amcset.exe [--seeds 8] [--scale 1] [--threads 2]

For each quantity it prints the means of both modes, their relative
difference and z = (fast - exact) / sqrt(s_e^2/n + s_f^2/n), where s is the
seed-to-seed standard deviation (independent seeds: exact uses seeds
1..n, fast n+1..2n). With the same seeds the two modes share most of their
histories, so the paired relative difference (same seed, fast vs exact)
shows the systematic effect of the approximations directly. |z| > 3 is
expected for about 0.3 % of the quantities by chance. Standard library only.
"""

import argparse
import json
import math
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench import CASES, DEFAULT_CASES  # noqa: E402


def leaves(o, p=""):
    if isinstance(o, dict):
        for k, v in o.items():
            yield from leaves(v, f"{p}.{k}" if p else k)
    elif isinstance(o, list):
        for i, v in enumerate(o):
            yield from leaves(v, f"{p}[{i}]")
    elif isinstance(o, (int, float)) and not isinstance(o, bool):
        yield p, float(o)


SKIP = ("elapsed", "seconds", "seed", "ions", "energy_keV", "mass_amu", "angle_deg", ".z", "first_ion")


def run(exe, case_dir, ions, seed, fast, out, threads):
    cmd = [exe, "trim", str(case_dir / "TRIM.IN"), "-n", str(ions), "-s", str(seed), "-j", str(threads), "-q",
           "--formats", "trim", "-o", str(out)] + (["--fast"] if fast else [])
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{case_dir.name} seed {seed}: exit {p.returncode}\n{p.stderr}")
    s = json.loads((out / "summary.json").read_text())
    return {k: v for k, v in leaves(s) if not any(x in k for x in SKIP)}


def mean_sd(v):
    m = sum(v) / len(v)
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else 0.0
    return m, sd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--cases", default=str(DEFAULT_CASES))
    ap.add_argument("--seeds", type=int, default=8)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply every ion count")
    ap.add_argument("--only", nargs="*", help="substrings of the case names to run")
    ap.add_argument("--threads", type=int, default=2, help="threads per run (results do not depend on it)")
    ap.add_argument("--json", help="write every comparison to this file")
    args = ap.parse_args()
    n = args.seeds
    rows = []
    with tempfile.TemporaryDirectory() as tmp:
        for name, ions in CASES:
            if args.only and not any(o in name for o in args.only):
                continue
            case_dir = Path(args.cases) / name
            ions = max(1, int(ions * args.scale))
            res = {"exact": [], "fast": [], "paired": []}
            for i in range(n):
                e = run(args.exe, case_dir, ions, 1 + i, False, Path(tmp) / f"{name}_e{i}", args.threads)
                f = run(args.exe, case_dir, ions, 1 + n + i, True, Path(tmp) / f"{name}_f{i}", args.threads)
                fp = run(args.exe, case_dir, ions, 1 + i, True, Path(tmp) / f"{name}_p{i}", args.threads)
                res["exact"].append(e)
                res["fast"].append(f)
                res["paired"].append((e, fp))
            for k in res["exact"][0]:
                if not all(k in r for r in res["exact"] + res["fast"]):
                    continue
                me, se = mean_sd([r[k] for r in res["exact"]])
                mf, sf = mean_sd([r[k] for r in res["fast"]])
                err = math.sqrt((se * se + sf * sf) / n)
                z = (mf - me) / err if err > 0 else 0.0
                pe = sum(a[k] for a, _ in res["paired"])
                pf = sum(b[k] for _, b in res["paired"])
                paired = (pf - pe) / abs(pe) if pe != 0 else 0.0
                rows.append({"case": name, "ions": ions, "quantity": k, "exact": me, "fast": mf,
                             "rel": (mf - me) / abs(me) if me != 0 else 0.0, "z": z, "paired_rel": paired})
            worst = max((r for r in rows if r["case"] == name), key=lambda r: abs(r["z"]), default=None)
            if worst:
                print(f"{name:40s} n={ions:6d} x {n} seeds: worst |z| {abs(worst['z']):.2f} ({worst['quantity']})",
                      flush=True)
    zs = [abs(r["z"]) for r in rows if r["z"] != 0.0]
    big = [r for r in rows if abs(r["z"]) > 3]
    paired = sorted((abs(r["paired_rel"]), r["case"], r["quantity"]) for r in rows)
    print(f"\n{len(rows)} quantities in {len({r['case'] for r in rows})} cases; |z| > 2: "
          f"{sum(1 for z in zs if z > 2)}, |z| > 3: {len(big)} (by chance about {0.0027 * len(zs):.1f})")
    if zs:
        print(f"median |z| {sorted(zs)[len(zs) // 2]:.2f} (0.67 expected)")
    for r in big:
        print(f"  |z| {abs(r['z']):.2f}  {r['case']}  {r['quantity']}: exact {r['exact']:.6g} fast {r['fast']:.6g}")
    print("largest paired relative differences (same seeds):")
    for p, c, q in paired[-8:][::-1]:
        print(f"  {p:.2e}  {c}  {q}")
    if args.json:
        Path(args.json).write_text(json.dumps(rows, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())

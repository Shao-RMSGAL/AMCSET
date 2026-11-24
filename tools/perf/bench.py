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

"""Speed benchmark and bit-identity check of AMCSET's SRIM mode.

Runs a fixed set of TRIM reference inputs (TRIM.IN of the validation cases)
with a given amcset executable and records the calculation time that
summary.json reports (``elapsed_seconds``: the transport only, without
start-up and output) and the wall time of the whole process.

    python tools/perf/bench.py run --exe build/Release/amcset.exe --tag after
    python tools/perf/bench.py run --exe build/perf/baseline/amcset.exe --tag before
    python tools/perf/bench.py compare build/perf/runs/before build/perf/runs/after
    python tools/perf/bench.py table build/perf/runs/before build/perf/runs/after
    python tools/perf/bench.py ab --exe old/amcset.exe --tag old --exe new/amcset.exe --tag new

``compare`` checks that two runs wrote identical output files (lines with
time stamps, directories and run times are ignored), which is how a pure
optimisation is shown not to change any result. ``table`` prints the speed-up
per case. Only the Python standard library is used.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

DEFAULT_CASES = Path(os.environ.get("USERPROFILE", str(Path.home()))) / \
    "Applications" / "SRIM-2013-Pro" / "_amcset" / "results" / "trim3"

# (case folder, ions for the single-thread run). Sized for a few seconds per
# case on one core of a Ryzen 7 7800X3D with the 0.1.0 code.
CASES = [
    ("434_rep_fe_fe_1keV_fc_s1", 25000),
    ("422_rep_fe_fe_10keV_fc_s1", 6000),
    ("377_rep_fe_fe_100keV_fc_s1", 1000),
    ("328_rep_fe_fe_1000keV_fc_s1", 160),
    ("330_rep_fe_316ss_2000keV_fc_s1", 120),
    ("324_rep_fe_fe_1000keV_kp_ed25eb3_s1", 12000),
    ("341_rep_h_fe_1000keV_kp_s1", 100000),
    ("321_opt_he_w_2000keV_kp", 36000),
    ("426_rep_ar_cu_1keV_sputt_s1", 25000),
    ("355_lay_fe_crfe_order_2000keV_fc_s1", 120),
    ("385_gas_ar_ar_100keV_fc_gas", 12),
    ("386_dat_fe_fe_ions_fc", 400),
    ("375_dat_fe_fe_pka_fc", 3000),
]

# Lines that legitimately differ between runs (time stamps, output directory,
# run times) or between builds (the scenario key physics.fast_math, added
# after 0.1.0, appears in results.json only in newer builds).
VOLATILE = ("created", "\"directory\"", "elapsed", "Calculation time", "\"fast_math\": false")


def run_case(exe, case_dir, out_dir, ions, seed, threads, formats, extra=""):
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(exe), "trim", str(case_dir / "TRIM.IN"), "-n", str(ions), "-s", str(seed),
           "-j", str(threads), "-q", "--formats", formats, "-o", str(out_dir)] + extra.split()
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    if proc.returncode != 0:
        raise RuntimeError(f"{case_dir.name}: exit {proc.returncode}\n{proc.stderr}")
    summary = json.loads((out_dir / "summary.json").read_text())
    return {"case": case_dir.name, "ions": ions, "seed": seed, "threads": threads,
            "elapsed": summary["elapsed_seconds"], "wall": wall}


def cmd_run(args):
    cases_root = Path(args.cases)
    out_root = Path(args.out) / args.tag
    selected = [c for c in CASES if not args.only or any(s in c[0] for s in args.only)]
    results = []
    for name, ions in selected:
        n = max(1, int(ions * args.scale))
        if args.threads != 1:
            n = max(n, 16 * args.threads)  # enough blocks for every thread
        for seed in args.seeds:
            best = None
            for rep in range(args.repeat):
                r = run_case(args.exe, cases_root / name, out_root / name / str(seed), n, seed,
                             args.threads, args.formats, args.args)
                best = r if best is None or r["elapsed"] < best["elapsed"] else best
            results.append(best)
            print(f"{name:40s} n={n:7d} seed={seed:4d}  calc {best['elapsed']:8.3f} s"
                  f"  wall {best['wall']:8.3f} s", flush=True)
    (out_root / "timings.json").write_text(json.dumps(results, indent=1))


def cmd_ab(args):
    """Interleaved A/B timing: for each case and repetition run every
    executable in turn (drift in clock speed or background load then hits
    all of them alike) and keep each one's fastest run."""
    cases_root = Path(args.cases)
    exes = [Path(e) for e in args.exe]
    tags = args.tag or [f"exe{i}" for i in range(len(exes))]
    selected = [c for c in CASES if not args.only or any(s in c[0] for s in args.only)]
    results = {t: [] for t in tags}
    for name, ions in selected:
        n = max(1, int(ions * args.scale))
        if args.threads != 1:
            n = max(n, 16 * args.threads)
        best = {}
        for rep in range(args.repeat):
            for exe, tag, extra in zip(exes, tags, args.args or [""] * len(exes)):
                r = run_case(exe, cases_root / name, Path(args.out) / tag / name / str(args.seed), n,
                             args.seed, args.threads, args.formats, extra)
                if tag not in best or r["elapsed"] < best[tag]["elapsed"]:
                    best[tag] = r
        line = f"{name:40s} n={n:7d}"
        for tag in tags:
            results[tag].append(best[tag])
            line += f"  {tag} {best[tag]['elapsed']:7.3f}/{best[tag]['wall']:7.3f}"
        if len(tags) > 1:
            line += f"  x{best[tags[0]]['elapsed'] / max(best[tags[-1]]['elapsed'], 1e-9):5.2f}"
            line += f" (wall x{best[tags[0]]['wall'] / max(best[tags[-1]]['wall'], 1e-9):5.2f})"
        print(line, flush=True)
    for tag in tags:
        d = Path(args.out) / tag
        d.mkdir(parents=True, exist_ok=True)
        (d / "timings.json").write_text(json.dumps(results[tag], indent=1))


def stable_lines(path):
    with open(path, "rb") as f:
        data = f.read()
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError:
        return [data]
    # Paths may be spelled with either separator (the data directory is found
    # relative to the executable).
    return [ln.replace("\\", "/") for ln in text.splitlines() if not any(v in ln for v in VOLATILE)]


def cmd_compare(args):
    a, b = Path(args.a), Path(args.b)
    files = sorted(p.relative_to(a) for p in a.rglob("*") if p.is_file() and p.name != "timings.json")
    same = differ = missing = 0
    for rel in files:
        pb = b / rel
        if not pb.exists():
            missing += 1
            print(f"missing in {b}: {rel}")
            continue
        if stable_lines(a / rel) == stable_lines(pb):
            same += 1
        else:
            differ += 1
            print(f"DIFFERS: {rel}")
    print(f"{same} identical, {differ} different, {missing} missing")
    return 1 if differ or missing else 0


def cmd_table(args):
    runs = [json.loads((Path(d) / "timings.json").read_text()) for d in args.dirs]
    base = {(r["case"], r["seed"]): r for r in runs[0]}
    head = f"{'case':40s} {'ions':>7s} " + " ".join(f"{Path(d).name:>12s}" for d in args.dirs)
    if len(runs) > 1:
        head += "  speed-up"
    print(head)
    logsum, count = 0.0, 0
    for key, r0 in base.items():
        line = f"{key[0]:40s} {r0['ions']:7d} {r0['elapsed']:12.3f}"
        last = r0
        for run in runs[1:]:
            other = {(r["case"], r["seed"]): r for r in run}.get(key)
            line += f" {other['elapsed']:12.3f}" if other else f" {'-':>12s}"
            last = other or last
        if len(runs) > 1 and last is not r0 and last["elapsed"] > 0:
            s = r0["elapsed"] / last["elapsed"]
            import math
            logsum += math.log(s)
            count += 1
            line += f"  {s:7.2f}x"
        print(line)
    if count:
        import math
        print(f"geometric mean speed-up: {math.exp(logsum / count):.2f}x over {count} cases")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--exe", required=True)
    r.add_argument("--tag", required=True)
    r.add_argument("--cases", default=str(DEFAULT_CASES))
    r.add_argument("--out", default="build/perf/runs")
    r.add_argument("--threads", type=int, default=1)
    r.add_argument("--seeds", type=int, nargs="+", default=[101])
    r.add_argument("--repeat", type=int, default=1, help="keep the fastest of N runs")
    r.add_argument("--scale", type=float, default=1.0, help="multiply every ion count")
    r.add_argument("--only", nargs="*", help="substrings of the case names to run")
    r.add_argument("--formats", default="trim,csv,json")
    r.add_argument("--args", default="", help="extra amcset options, e.g. --args=--fast")
    r.set_defaults(func=cmd_run)
    ab = sub.add_parser("ab")
    ab.add_argument("--exe", action="append", required=True)
    ab.add_argument("--tag", action="append")
    ab.add_argument("--cases", default=str(DEFAULT_CASES))
    ab.add_argument("--out", default="build/perf/runs")
    ab.add_argument("--threads", type=int, default=1)
    ab.add_argument("--seed", type=int, default=101)
    ab.add_argument("--repeat", type=int, default=3)
    ab.add_argument("--scale", type=float, default=1.0)
    ab.add_argument("--only", nargs="*")
    ab.add_argument("--formats", default="trim,csv,json")
    ab.add_argument("--args", action="append",
                    help="extra amcset options per --exe, in order (e.g. --args= --args=--fast)")
    ab.set_defaults(func=cmd_ab)
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    c.set_defaults(func=cmd_compare)
    t = sub.add_parser("table")
    t.add_argument("dirs", nargs="+")
    t.set_defaults(func=cmd_table)
    args = p.parse_args()
    sys.exit(args.func(args) or 0)


if __name__ == "__main__":
    main()

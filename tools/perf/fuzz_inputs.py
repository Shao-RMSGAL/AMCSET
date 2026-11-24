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

"""Robustness check of the input readers: feed amcset truncated and mutated
TRIM.IN files and scenario JSON files and look for crashes and hangs.

    python tools/perf/fuzz_inputs.py --exe build/Release/amcset.exe --out build/fuzz [--count 400]

Every input goes through ``amcset validate``; inputs that validate are also
run with 3 ions. A crash (exit code 3, with a crash report) or a time-out is a
bug; exit codes 1 and 2 (an error message) are the expected outcome for bad
input. Uses only the standard library.
"""

import argparse
import json
import random
import subprocess
import sys
from pathlib import Path

NUMBERS = ["0", "-1", "-0", "1e308", "-1e308", "1e-320", "nan", "inf", "99999999999999999999",
           "4294967296", "0.0000001", "", "x", "1,5", "3.2.1"]


def mutate_text(text, rng):
    lines = text.splitlines()
    kind = rng.randrange(6)
    if kind == 0 and lines:  # truncate
        return "\n".join(lines[: rng.randrange(len(lines))]) + "\n"
    if kind == 1 and lines:  # drop a line
        i = rng.randrange(len(lines))
        return "\n".join(lines[:i] + lines[i + 1:]) + "\n"
    if kind == 2 and lines:  # duplicate a line
        i = rng.randrange(len(lines))
        return "\n".join(lines[: i + 1] + lines[i:]) + "\n"
    if kind == 3:  # replace one token with a hostile number
        toks = text.split(" ")
        idx = [i for i, t in enumerate(toks) if any(c.isdigit() for c in t)]
        if idx:
            toks[rng.choice(idx)] = rng.choice(NUMBERS)
        return " ".join(toks)
    if kind == 4 and text:  # flip a byte
        i = rng.randrange(len(text))
        return text[:i] + chr(rng.randrange(32, 127)) + text[i + 1:]
    return text * 2  # concatenated copies


def mutate_json(obj, rng):
    """Replace one leaf with a hostile value."""
    paths = []

    def walk(o, p):
        if isinstance(o, dict):
            for k, v in o.items():
                walk(v, p + [k])
        elif isinstance(o, list):
            for i, v in enumerate(o):
                walk(v, p + [i])
        else:
            paths.append(p)

    walk(obj, [])
    if not paths:
        return obj
    p = rng.choice(paths)
    target = obj
    for k in p[:-1]:
        target = target[k]
    old = target[p[-1]]
    hostile = [0, -1, 1e308, -1e308, 2**64, "", "0 eV", "-5 eV", "1e308 MeV", "1 kg", "nan",
               [], {}, None, True, "Xx"]
    if isinstance(old, str) and " " in old:
        num, unit = old.split(" ", 1)
        hostile += [f"{n} {unit}" for n in NUMBERS]
    target[p[-1]] = rng.choice(hostile)
    return obj


def run(cmd, timeout):
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        return p.returncode, p.stderr[-2000:]
    except subprocess.TimeoutExpired:
        return "timeout", ""


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True)
    ap.add_argument("--out", default="build/fuzz")
    ap.add_argument("--count", type=int, default=400)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--trim-in", nargs="*", default=[], help="TRIM.IN files to start from")
    args = ap.parse_args()
    rng = random.Random(args.seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    template = subprocess.run([args.exe, "template"], capture_output=True, text=True).stdout
    base_json = json.loads(template)
    trim_ins = [Path(p).read_text(errors="replace") for p in args.trim_in]
    problems = []
    gaps = []
    counts = {}
    for i in range(args.count):
        if trim_ins and i % 2 == 0:
            name = out / f"case{i:04d}" / "TRIM.IN"
            name.parent.mkdir(exist_ok=True)
            name.write_text(mutate_text(rng.choice(trim_ins), rng))
        else:
            name = out / f"case{i:04d}.json"
            obj = json.loads(json.dumps(base_json))
            obj["ions"] = 3
            if rng.random() < 0.25:  # the automatic target depth (quick test runs)
                obj["auto_depth"] = {"enabled": True, "test_ions": 20}
            for _ in range(rng.randrange(1, 3)):
                obj = mutate_json(obj, rng)
            text = json.dumps(obj)
            if rng.random() < 0.2:
                text = mutate_text(text, rng)
            name.write_text(text)
        code, err = run([args.exe, "validate", str(name)], 20)
        counts[f"validate:{code}"] = counts.get(f"validate:{code}", 0) + 1
        if code == 3 or code == "timeout" or (isinstance(code, int) and code not in (0, 1, 2)):
            problems.append((str(name), "validate", code, err))
            continue
        if code == 0:
            # (A TRIM.IN has no auto_depth key; --auto-depth sizes some of them.)
            extra = ["--auto-depth"] if name.suffix != ".json" and rng.random() < 0.25 else []
            cmd = [args.exe, "run" if name.suffix == ".json" else "trim", str(name), *extra, "-n", "3", "-j", "2",
                   "-q", "--formats", "trim", "-o", str(name.parent / (name.stem + "_out"))]
            code, err = run(cmd, 120)
            counts[f"run:{code}"] = counts.get(f"run:{code}", 0) + 1
            if code == 3 or code == "timeout" or (isinstance(code, int) and code not in (0, 1, 2)):
                problems.append((str(name), "run", code, err))
            elif code != 0:
                # Valid by `validate` but refused by `run`: a gap in the
                # validation unless the cause is outside the file (a missing
                # stopping table or SRIM directory).
                gaps.append((str(name), err.strip().splitlines()[-1] if err.strip() else ""))
    print("outcomes:", dict(sorted(counts.items(), key=lambda kv: str(kv[0]))))
    for g in gaps:
        print("VALID BUT NOT RUN:", g[0], "-", g[1])
    for p in problems:
        print("PROBLEM:", p[0], p[1], p[2])
        print(p[3])
    print(f"{len(problems)} problems in {args.count} inputs")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())

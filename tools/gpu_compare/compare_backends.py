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

"""AMCSET on the CPU against AMCSET on GPU backends, statistically.

A GPU runs the CPU's trim2013 transport, but its math library differs in the
last bit (and mixed precision more), so the histories part: the results can
only agree within their statistics. This script runs a set of cases on the
CPU and on each backend with independent seeds, and for every quantity of
summary.json computes z = (mean_gpu - mean_cpu) / sqrt(se_gpu^2 + se_cpu^2)
over the runs. With the same physics and n runs on each side, z follows a
t distribution with about 2(n - 1) degrees of freedom, as the standard errors
are themselves estimated: for 6 runs about 7.3 % beyond |z| = 2 and 1.3 %
beyond 3 (a normal distribution: 4.6 % and 0.27 %). The quantities of one
case share their runs and the CPU's runs are shared by every backend, so
deviations come in correlated groups; check a group with more and larger
runs before calling it real. docs/GPU.md "CPU against GPU".

  python tools/gpu_compare/compare_backends.py --amcset build/Release/amcset.exe \
      --backends opencl:double,opencl:mixed,hip:mixed --runs 6 --out build/gpu_compare
"""

import argparse
import json
import math
import os
import subprocess
import sys
import time

# The cases: (name, scenario, ions per run). Built-in physics: the SRIM-2013
# tables that come with AMCSET.
def fe(energy, damage, thickness):
    return {"ion": {"element": "Fe", "energy": energy},
            "layers": [{"name": "Fe", "thickness": thickness, "density": "7.8658 g/cm3",
                        "elements": [{"element": "Fe", "Ed": "40 eV", "Eb": "0 eV"}]}],
            "damage": damage}

CASES = [
    ("fe50keV_full", fe("50 keV", "full", "200 nm"), 20000),
    ("fe50keV_kp", fe("50 keV", "kp", "200 nm"), 100000),
    ("fe1MeV_full", fe("1 MeV", "full", "1500 nm"), 2000),
    ("ar30keV_SiO2_Si_full", {
        "ion": {"element": "Ar", "energy": "30 keV", "angle": "7 deg"},
        "layers": [{"name": "SiO2", "thickness": "15 nm", "density": "2.2 g/cm3",
                    "elements": [{"element": "Si", "stoichiometry": 1}, {"element": "O", "stoichiometry": 2}]},
                   {"name": "Si", "thickness": "200 nm", "density": "2.3212 g/cm3", "elements": [{"element": "Si"}]}],
        "damage": "full"}, 20000),
    ("h100keV_Si_full", {
        "ion": {"element": "H", "energy": "100 keV"},
        "layers": [{"name": "Si", "thickness": "1500 nm", "density": "2.3212 g/cm3", "elements": [{"element": "Si"}]}],
        "damage": "full"}, 50000),
    ("au200keV_Si_full", {
        "ion": {"element": "Au", "energy": "200 keV"},
        "layers": [{"name": "Si", "thickness": "300 nm", "density": "2.3212 g/cm3", "elements": [{"element": "Si"}]}],
        "damage": "full"}, 4000),
    ("ar2keV_Cu_monolayer", {
        "ion": {"element": "Ar", "energy": "2 keV", "angle": "45 deg"},
        "layers": [{"name": "Cu", "thickness": "20 nm", "density": "8.96 g/cm3", "elements": [{"element": "Cu"}]}],
        "damage": "monolayer"}, 20000),
]

SECTIONS = ["range", "backscattered", "transmitted", "per_ion", "energy_percent", "damage_models"]
SKIP = {"count", "stopped_ions", "dpa_fluence_per_cm2"}


def quantities(summary):
    """The numbers of summary.json to compare, by name."""
    q = {}
    for sec in SECTIONS:
        for k, v in summary.get(sec, {}).items():
            if k not in SKIP and isinstance(v, (int, float)):
                q[f"{sec}.{k}"] = float(v)
    for el in summary.get("elements", []):
        name = el.get("name", el.get("symbol", "?"))
        for k, v in el.items():
            if isinstance(v, (int, float)) and k not in ("z", "mass_amu", "index"):
                q[f"element.{name}.{k}"] = float(v)
    return q


def run(amcset, scenario, ions, seed, backend, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "scenario.json")
    sc = dict(scenario)
    sc["stopping"] = {"model": "srim"}
    with open(path, "w") as f:
        json.dump(sc, f)
    args = [amcset, "run", path, "-q", "-n", str(ions), "-s", str(seed), "-o", out_dir, "--formats", "trim",
            "--set", "output.exit_lists=false", "--set", "output.range_3d=false"]
    if backend != "cpu":
        name, precision = backend.split(":")
        args += ["--compute", name, "--precision", precision]
    t0 = time.time()
    p = subprocess.run(args, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(f"{' '.join(args)}\n{p.stdout}\n{p.stderr}")
    with open(os.path.join(out_dir, "summary.json")) as f:
        return quantities(json.load(f)), time.time() - t0


def t_tail(x, dof):
    """P(|t| > x) for Student's t with dof degrees of freedom (Simpson's rule)."""
    def pdf(t):
        return math.exp(math.lgamma((dof + 1) / 2) - math.lgamma(dof / 2)
                        - (dof + 1) / 2 * math.log1p(t * t / dof)) / math.sqrt(dof * math.pi)
    n = 2000
    h = x / n
    inner = pdf(0) + pdf(x) + sum((4 if i % 2 else 2) * pdf(i * h) for i in range(1, n))
    return 1.0 - 2.0 * inner * h / 3.0


def summarize(all_z, runs):
    """Per backend: how the z values compare with Student's t of 2(runs - 1) degrees of freedom."""
    dof = 2 * (runs - 1)
    expect2 = t_tail(2.0, dof) if dof > 0 else float("nan")
    expect3 = t_tail(3.0, dof) if dof > 0 else float("nan")
    rms_expected = math.sqrt(dof / (dof - 2)) if dof > 2 else float("inf")
    summary = {}
    for b, zs in all_z.items():
        n = len(zs)
        if n == 0:
            continue
        over2 = sum(1 for z in zs if abs(z) > 2) / n
        over3 = sum(1 for z in zs if abs(z) > 3)
        mean = sum(zs) / n
        rms = math.sqrt(sum(z * z for z in zs) / n)
        print(f"{b:15s} {n} quantities: mean z {mean:+.2f}, rms z {rms:.2f} ({rms_expected:.2f} expected), "
              f"{100 * over2:.1f} % beyond |z| = 2 ({100 * expect2:.1f} % expected), "
              f"{over3} beyond 3 ({n * expect3:.1f} expected)")
        summary[b] = {"n": n, "mean_z": mean, "rms_z": rms, "frac_over_2": over2, "n_over_3": over3, "dof": dof,
                      "expected_frac_over_2": expect2, "expected_n_over_3": n * expect3,
                      "expected_rms_z": rms_expected}
    return summary


def resummarize(path):
    """Recompute the summary of an existing analysis.json (--resummarize)."""
    with open(path) as f:
        results = json.load(f)
    all_z = {}
    runs = 6
    for name, case in results.items():
        if name.startswith("_"):
            continue
        runs = case.get("runs", runs)
        for q in case["quantities"].values():
            for b, v in q.items():
                if b != "cpu" and math.isfinite(v[2]) and math.hypot(v[1], q["cpu"][1]) > 0.0:
                    all_z.setdefault(b, []).append(v[2])
    results["_summary"] = summarize(all_z, runs)
    with open(path, "w") as f:
        json.dump(results, f, indent=1)


def stats(values):
    n = len(values)
    m = sum(values) / n
    var = sum((v - m) ** 2 for v in values) / (n - 1) if n > 1 else 0.0
    return m, math.sqrt(var / n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--amcset", required=True)
    ap.add_argument("--backends", default="opencl:mixed")
    ap.add_argument("--runs", type=int, default=6)
    ap.add_argument("--scale", type=float, default=1.0, help="ions per run times this")
    ap.add_argument("--cases", default="", help="comma-separated case names (default: all)")
    ap.add_argument("--out", default="build/gpu_compare")
    ap.add_argument("--resummarize", metavar="ANALYSIS_JSON",
                    help="only recompute the summary of an existing analysis.json")
    a = ap.parse_args()
    if a.resummarize:
        resummarize(a.resummarize)
        return 0
    backends = ["cpu"] + [b for b in a.backends.split(",") if b]
    cases = [c for c in CASES if not a.cases or c[0] in a.cases.split(",")]
    results = {}
    lines = []
    all_z = {b: [] for b in backends[1:]}
    for name, scenario, ions in cases:
        ions = max(16, int(ions * a.scale))
        per_backend = {}
        for bi, b in enumerate(backends):
            runs, seconds = [], []
            for r in range(a.runs):
                seed = 1000 * (bi + 1) + r + 1  # independent seeds for every backend and run
                q, dt = run(a.amcset, scenario, ions, seed, b, os.path.join(a.out, name, b.replace(":", "_"), str(r)))
                runs.append(q)
                seconds.append(dt)
            per_backend[b] = {"runs": runs, "seconds": sum(seconds) / len(seconds)}
            print(f"{name:24s} {b:15s} {ions:7d} ions x {a.runs}: {per_backend[b]['seconds']:.2f} s per run", flush=True)
        cpu = per_backend["cpu"]["runs"]
        keys = sorted(set.intersection(*[set(r.keys()) for r in cpu]))
        results[name] = {"ions": ions, "runs": a.runs, "seconds": {b: per_backend[b]["seconds"] for b in backends},
                         "quantities": {}}
        for b in backends[1:]:
            gpu = per_backend[b]["runs"]
            worst = (0.0, "")
            for k in keys:
                mc, sc = stats([r[k] for r in cpu])
                mg, sg = stats([r[k] for r in gpu if k in r])
                se = math.sqrt(sc * sc + sg * sg)
                if se == 0.0:
                    z = 0.0 if mc == mg else float("inf")
                else:
                    z = (mg - mc) / se
                results[name]["quantities"].setdefault(k, {"cpu": [mc, sc]})[b] = [mg, sg, z]
                if math.isfinite(z) and se > 0.0:  # a constant (e.g. no transmission) says nothing
                    all_z[b].append(z)
                if abs(z) > abs(worst[0]):
                    worst = (z, k)
            lines.append(f"{name:24s} {b:15s} max |z| {abs(worst[0]):5.2f} ({worst[1]})")
    print()
    print("Largest deviation per case and backend:")
    for line in lines:
        print("  " + line)
    print()
    results["_summary"] = summarize(all_z, a.runs)
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "analysis.json"), "w") as f:
        json.dump(results, f, indent=1)
    print(f"\nWritten: {os.path.join(a.out, 'analysis.json')}")


if __name__ == "__main__":
    sys.exit(main())

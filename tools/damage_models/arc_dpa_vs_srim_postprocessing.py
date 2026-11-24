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

"""AMCSET's arc-dpa against published SRIM-2013 post-processing.

Mitsi, Koutsomitis and Apostolopoulos, "On the use of SRIM for calculating
arc-dpa exposure", Nucl. Instrum. Methods B (2023), doi
10.1016/j.nimb.2023.165145 (arXiv 2307.12867), applied the arc-dpa
efficiency of Nordlund et al. (2018) to every primary of SRIM-2013 quick
calculations (COLLISON.txt, their method M1; 10,000 ions). Their Table 5
gives, per primary, the NRT vacancies <nu_d> and <nu_d,arc>; Table 4 gives
per-ion totals for Fe into Fe. This script runs the same cases in AMCSET's
quick calculation with their Ed values and compares arc/NRT (the number of
primaries cancels) and the Fe per-ion totals.

    python tools/damage_models/arc_dpa_vs_srim_postprocessing.py --exe build/Release/amcset.exe \
        --out build/damage_models/mitsi2023 [--ions 2000] [--seeds 2]
"""

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

# Their Table 2: Ed (eV); arc constants b, c of Nordlund et al. 2018. Density g/cm3.
TARGETS = {
    "Fe": (40.0, -0.568, 0.286, 7.8658),
    "Ni": (40.0, -1.01, 0.23, 8.90),
    "Cu": (29.0, -0.68, 0.16, 8.92),
    "Pd": (41.0, -0.88, 0.15, 12.0),
    "W": (90.0, -0.56, 0.12, 19.35),
    "Pt": (44.0, -1.12, 0.11, 21.45),
}
BEAMS = [("H", 1000.0), ("He", 1000.0), ("Al", 3000.0), ("Fe", 5000.0), ("Au", 10000.0)]
# Table 5, method M1: (<nu_d>, <nu_d,arc>) per primary.
TABLE5 = {
    "Fe": [(1.91, 1.38), (2.47, 1.62), (10.4, 4.24), (17.9, 6.56), (44.8, 14.8)],
    "Ni": [(1.92, 1.21), (2.54, 1.35), (10.6, 3.21), (18.6, 5.04), (47.7, 11.7)],
    "Cu": [(2.02, 1.27), (2.76, 1.45), (12.6, 3.29), (22.3, 4.95), (60.4, 11.4)],
    "Pd": [(1.89, 1.17), (2.47, 1.27), (11.6, 2.70), (22.7, 4.39), (67.0, 11.1)],
    "W": [(1.61, 1.20), (1.94, 1.29), (8.06, 2.51), (16.7, 3.85), (52.8, 8.88)],
    "Pt": [(1.78, 1.06), (2.37, 1.11), (11.7, 2.08), (25.6, 3.59), (87.4, 10.4)],
}
# Table 4, Fe into Fe, method M1: (keV, N_d, N_d,arc) per ion.
TABLE4 = [(78.7, 530.0, 217.0), (5000.0, 7870.0, 2900.0)]


def run(exe, work, target, ion, energy_keV, ions, seed):
    ed, b, c, density = TARGETS[target]
    d = work / f"{ion}_{energy_keV:g}keV_{target}_s{seed}"
    summary = d / "out" / "summary.json"
    if not summary.exists():
        d.mkdir(parents=True, exist_ok=True)
        sc = {
            "mode": "srim",
            "title": f"{ion} {energy_keV} keV into {target}",
            "ion": {"element": ion, "energy": f"{energy_keV} keV"},
            # Thick enough to stop every ion (1 MeV H in Fe: about 6 um).
            "layers": [{"name": target, "thickness": "40 um", "density": f"{density} g/cm3",
                        "elements": [{"element": target, "Ed": f"{ed} eV", "Eb": "0 eV",
                                      "damage_efficiency": {"arc_b": b, "arc_c": c, "rpa_b": "1000 eV",
                                                            "rpa_c": 1.0}}]}],
            "ions": ions,
            "seed": seed,
            "damage": "kp",
            "output": {"directory": str((d / "out").resolve()), "range_3d": False, "exit_lists": False,
                       "files": {"formats": ["json"]}},
        }
        (d / "scenario.json").write_text(json.dumps(sc), encoding="utf-8")
        r = subprocess.run([str(exe), "run", str(d / "scenario.json")], capture_output=True, text=True)
        if r.returncode != 0 or not summary.exists():
            sys.exit(f"{d.name}: amcset failed\n{r.stdout}\n{r.stderr}")
    s = json.loads(summary.read_text(encoding="utf-8"))
    dm = s["damage_models"]
    return dm["nrt_per_ion"], dm["arc_dpa_per_ion"], s["per_ion"]["vacancies"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--ions", type=int, default=2000)
    ap.add_argument("--seeds", type=int, default=2)
    a = ap.parse_args()
    work = a.out / "runs"
    lines = ["arc/NRT over all primaries: AMCSET quick calculation vs Mitsi et al. 2023 Table 5 (M1)",
             f"{a.ions} ions x {a.seeds} seeds per case; +- standard error over seeds", "",
             f"{'target':6s} {'beam':12s} {'AMCSET arc/NRT':>18s} {'Mitsi':>7s} {'diff %':>7s}"]
    diffs = []
    for target in TARGETS:
        for (ion, e), (nu, nu_arc) in zip(BEAMS, TABLE5[target]):
            ratios = []
            for k in range(a.seeds):
                nrt, arc, vac = run(a.exe, work, target, ion, e, a.ions, 101 + k)
                if abs(nrt - vac) > 1e-6 * vac:
                    sys.exit(f"{target} {ion}: NRT {nrt} differs from TRIM's quick vacancies {vac}")
                ratios.append(arc / nrt)
            m = sum(ratios) / len(ratios)
            se = (math.sqrt(sum((x - m) ** 2 for x in ratios) / (len(ratios) - 1) / len(ratios))
                  if len(ratios) > 1 else float("nan"))
            ref = nu_arc / nu
            diffs.append(100.0 * (m / ref - 1.0))
            lines.append(f"{target:6s} {ion + ' ' + format(e / 1000, 'g') + ' MeV':12s} {m:10.4f} +-{se:6.4f} "
                         f"{ref:7.4f} {diffs[-1]:7.2f}")
    rms = math.sqrt(sum(d * d for d in diffs) / len(diffs))
    lines += ["", f"rms difference {rms:.2f} %, largest {max(diffs, key=abs):.2f} % (Table 5 values have 3 digits)",
              "", "Fe into Fe per ion: AMCSET vs Mitsi et al. 2023 Table 4 (M1)",
              f"{'keV':>7s} {'NRT':>9s} {'ref':>7s} {'arc':>9s} {'ref':>7s}"]
    for e, nd, nd_arc in TABLE4:
        res = [run(a.exe, work, "Fe", "Fe", e, a.ions, 101 + k) for k in range(a.seeds)]
        nrt = sum(r[0] for r in res) / len(res)
        arc = sum(r[1] for r in res) / len(res)
        lines.append(f"{e:7g} {nrt:9.1f} {nd:7.0f} {arc:9.1f} {nd_arc:7.0f}")
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "comparison.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

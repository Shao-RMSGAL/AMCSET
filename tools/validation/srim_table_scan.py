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

"""Scan the bundled SRIM-2013 stopping tables for ions out of line with their neighbours.

Electronic stopping at a given velocity changes smoothly with the ion's atomic
number Z1 except for the Z1 oscillations at low velocity (at most some 30 %
between neighbours, largest below ~0.1 MeV/u). An ion whose stopping departs
from the mean of its two neighbours by much more than that, over a wide band of
velocities, is suspect. The known case is Ni: SRIM-2013 overestimates the
ranges of Ni ions by 1.3-1.5 times (Gurbich et al., Appl. Radiat. Isot. 215,
111590 (2025)).

    python tools/validation/srim_table_scan.py --exe build/Release/amcset.exe --out build/validation/table_scan
        [--targets C,Si,Fe,Ni,Cu,Mo,Ag,W,Au] [--threshold 0.2]

For every target and every ion Z1 = 3..90 the ratio of S_e(Z1) to the
geometric mean of its neighbours' S_e at the same energy per mass (0.01-10
MeV/u, 10 points per decade), once with Z1 +- 1 and once with Z1 +- 2. The
geometric mean follows the smooth Z1 dependence (S_e grows like a power of
Z1); requiring both ratios to leave the band keeps an ion from being listed
only because a neighbour is out of line. Writes OUT/ratios.csv and
OUT/results.txt (the ions whose two ratios both leave 1 +- threshold at the
same velocity of 0.1 MeV/u or more, by target). Standard library only.
"""

import argparse
import math
import csv
import subprocess
import sys
from pathlib import Path

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
           "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
           "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()


def table(exe, ion, target):
    """{energy per amu (keV/u, rounded): S_e (eV/A)} of the bundled SRIM tables."""
    r = subprocess.run([str(exe), "stopping", "--ion", ion, "--target", target, "--from", "1 keV", "--to",
                        "10 GeV", "--per-decade", "10", "--csv"], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"amcset stopping --ion {ion} --target {target} failed:\n{r.stderr}")
    mass = None
    rows = {}
    for line in r.stdout.splitlines():
        if line.startswith("#"):
            # "# AMCSET stopping table: Ni (57.94 amu) in Fe, ..."
            if mass is None and "amu)" in line:
                mass = float(line.split("(")[1].split()[0])
            continue
        parts = line.split(",")
        try:
            e_keV, se = float(parts[0]), float(parts[1])
        except (ValueError, IndexError):
            continue
        rows[e_keV] = se
    if mass is None:
        sys.exit(f"no mass in the output for {ion}")
    return mass, rows


def at(rows, e_keV):
    """S_e at e_keV by log-log interpolation."""
    import math
    keys = sorted(rows)
    if e_keV <= keys[0] or e_keV >= keys[-1]:
        return float("nan")
    lo = max(k for k in keys if k <= e_keV)
    hi = min(k for k in keys if k >= e_keV)
    if lo == hi:
        return rows[lo]
    f = (math.log(e_keV) - math.log(lo)) / (math.log(hi) - math.log(lo))
    return math.exp(math.log(rows[lo]) + f * (math.log(rows[hi]) - math.log(rows[lo])))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--targets", default="C,Si,Fe,Ni,Cu,Mo,Ag,W,Au")
    ap.add_argument("--threshold", type=float, default=0.2)
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    speeds = [0.01 * 10 ** (k / 10) for k in range(31)]  # MeV/u, 0.01 .. 10
    out_rows = []
    lines = ["Bundled SRIM-2013 tables: S_e(Z1) / geometric mean of the neighbours' S_e at the same energy per mass",
             f"(Z1 +- 1 and Z1 +- 2). Listed: ions whose two ratios both leave 1 +- {a.threshold:g} on the same side",
             "at one velocity of 0.1-10 MeV/u (the largest such departure).", ""]
    for target in a.targets.split(","):
        tables = {}
        for sym in SYMBOLS[:92]:
            tables[sym] = table(a.exe, sym, target)
        lines.append(f"Target {target}:")
        flagged = 0
        for i in range(2, 90):
            sym = SYMBOLS[i]
            m, rows = tables[sym]
            worst = None
            for v in speeds:
                se = at(rows, v * m * 1000.0)
                r = []
                for d in (1, 2):
                    nb = [at(tables[SYMBOLS[j]][1], v * tables[SYMBOLS[j]][0] * 1000.0) for j in (i - d, i + d)]
                    r.append(se / math.sqrt(nb[0] * nb[1]))
                out_rows.append({"target": target, "ion": sym, "z1": i + 1, "MeV_per_u": f"{v:.4g}",
                                 "S_e_eV_per_A": f"{se:.5g}", "ratio_z1pm1": f"{r[0]:.4f}",
                                 "ratio_z1pm2": f"{r[1]:.4f}"})
                if v < 0.1 - 1e-9 or r[0] != r[0] or r[1] != r[1]:
                    continue
                both = all(abs(x - 1.0) > a.threshold for x in r) and (r[0] - 1.0) * (r[1] - 1.0) > 0
                if both and (worst is None or abs(r[0] - 1.0) > abs(worst[0] - 1.0)):
                    worst = (r[0], r[1], v)
            if worst is not None:
                flagged += 1
                lines.append(f"  {sym:2s} (Z1 {i + 1:2d})  ratio {worst[0]:.3f} (Z1+-1), {worst[1]:.3f} (Z1+-2) "
                             f"at {worst[2]:.3g} MeV/u")
        if not flagged:
            lines.append("  none")
        lines.append("")
    with (a.out / "ratios.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(out_rows[0]))
        w.writeheader()
        w.writerows(out_rows)
    (a.out / "results.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

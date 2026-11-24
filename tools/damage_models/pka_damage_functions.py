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

"""Damage functions of primary knock-on atoms (PKAs) by damage model.

Starts PKAs of one energy deep inside a thick pure target (a TRIM.DAT recoil
source, isotropic directions) and reads, per PKA, what AMCSET's summary.json
reports: TRIM's full-cascade vacancies and replacements, the damage energy
T_d (TRIM's Robinson partition) and the NRT, Kinchin-Pease, arc-dpa and rpa
values of the same primaries (docs/OUTPUT.md "dpa"). Several seeds give the
standard error of the cascade numbers.

    python tools/damage_models/pka_damage_functions.py --exe build/Release/amcset.exe \
        --out build/damage_models [--elements Fe,Cu] [--energies 0.1,1,10] [--pkas 500] [--seeds 4]

Writes OUT/damage_functions.csv (one row per element, Eb and energy) and
OUT/damage_functions.txt (a readable table). Only the standard library is
used. The arc-dpa / rpa constants and Ed are those of Nordlund et al.,
Nat. Commun. 9, 1084 (2018), Table 1, the MD fits arc-dpa was made from.
"""

import argparse
import csv
import json
import math
import os
import random
import subprocess
import sys
from pathlib import Path

# z, density (g/cm3), Ed (eV) of Nordlund et al. 2018 Table 1, Es (eV, TRIM's defaults).
ELEMENTS = {
    "Fe": (26, 7.8658, 40.0, 4.34),
    "Cu": (29, 8.92, 33.0, 3.52),
    "Ni": (28, 8.90, 39.0, 4.46),
    "Pd": (46, 12.0, 41.0, 3.91),
    "Pt": (78, 21.45, 42.0, 5.86),
    "W": (74, 19.35, 70.0, 8.68),
}
ENERGIES_KEV = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200]


def isotropic(rng):
    """A unit vector with an isotropic direction."""
    c = 2.0 * rng.random() - 1.0
    phi = 2.0 * math.pi * rng.random()
    s = math.sqrt(max(0.0, 1.0 - c * c))
    return [c, s * math.cos(phi), s * math.sin(phi)]


def scenario(symbol, energy_keV, pkas, seed, eb_eV, out_dir):
    z, density, ed, es = ELEMENTS[symbol]
    rng = random.Random(f"{symbol}-{energy_keV}-{seed}")
    # Thick target, PKAs in the middle: no surface or back-side losses up to 200 keV.
    thickness_um = 20.0
    depth_A = 0.5 * thickness_um * 1e4
    particles = [{"element": symbol, "energy": f"{energy_keV * 1000.0:.6g} eV", "x": f"{depth_A} A",
                  "direction": isotropic(rng)} for _ in range(pkas)]
    return {
        "mode": "srim",
        "title": f"{symbol} PKA {energy_keV} keV",
        "ion": {"element": symbol, "energy": f"{energy_keV} keV"},
        "layers": [{"name": symbol, "thickness": f"{thickness_um} um", "density": f"{density} g/cm3",
                    "elements": [{"element": symbol, "Ed": f"{ed} eV", "Eb": f"{eb_eV} eV", "Es": f"{es} eV"}]}],
        "source": {"type": "recoils", "particles": particles},
        "seed": seed,
        "threads": 0,
        "damage": "full",
        "output": {"directory": str(out_dir), "range_3d": False, "exit_lists": False,
                   "files": {"formats": ["json"]}},
    }


def run_case(exe, symbol, energy_keV, pkas, seed, eb_eV, work):
    d = work / f"{symbol}_{energy_keV:g}keV_eb{eb_eV:g}_s{seed}"
    d.mkdir(parents=True, exist_ok=True)
    summary = d / "out" / "summary.json"
    if not summary.exists():
        sc = d / "scenario.json"
        # The output directory is resolved against the scenario's directory: make it absolute.
        sc.write_text(json.dumps(scenario(symbol, energy_keV, pkas, seed, eb_eV, (d / "out").resolve())),
                      encoding="utf-8")
        r = subprocess.run([str(exe), "run", str(sc)], capture_output=True, text=True)
        if r.returncode != 0 or not summary.exists():
            sys.exit(f"{d.name}: amcset failed\n{r.stdout}\n{r.stderr}")
    s = json.loads(summary.read_text(encoding="utf-8"))
    dm = s["damage_models"]
    return {
        "td": dm["damage_energy_eV_per_ion"],
        "trim_vac": s["per_ion"]["vacancies"],
        "trim_repl": s["per_ion"]["replacements"],
        "nrt": dm["nrt_per_ion"],
        "kp": dm["kinchin_pease_per_ion"],
        "arc": dm.get("arc_dpa_per_ion", float("nan")),
        "rpa": dm.get("rpa_per_ion", float("nan")),
    }


def mean_sem(values):
    n = len(values)
    m = sum(values) / n
    if n < 2:
        return m, float("nan")
    var = sum((v - m) ** 2 for v in values) / (n - 1)
    return m, math.sqrt(var / n)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--elements", default=",".join(ELEMENTS))
    ap.add_argument("--energies", default=",".join(f"{e:g}" for e in ENERGIES_KEV), help="keV")
    ap.add_argument("--pkas", type=int, default=400, help="PKAs per run")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--eb", default="0", help="lattice binding energies to run, eV (comma list)")
    a = ap.parse_args()
    work = a.out / "runs"
    rows = []
    for symbol in a.elements.split(","):
        for eb in [float(x) for x in a.eb.split(",")]:
            for e in [float(x) for x in a.energies.split(",")]:
                # Fewer primaries at high energy (cascade time grows with T).
                n = max(50, int(a.pkas * min(1.0, 10.0 / e))) if e > 10 else a.pkas
                res = [run_case(a.exe, symbol, e, n, 101 + k, eb, work) for k in range(a.seeds)]
                row = {"element": symbol, "Ed_eV": ELEMENTS[symbol][2], "Eb_eV": eb, "E_keV": e,
                       "pkas": n * a.seeds}
                for key in res[0]:
                    m, se = mean_sem([r[key] for r in res])
                    row[key] = m
                    row[key + "_sem"] = se
                rows.append(row)
                print(f"{symbol:2s} Eb={eb:g} {e:7g} keV  Td={row['td']:9.1f}  TRIM-FC={row['trim_vac']:8.2f}"
                      f"  NRT={row['nrt']:8.2f}  arc={row['arc']:8.2f}", flush=True)
    a.out.mkdir(parents=True, exist_ok=True)
    with open(a.out / "damage_functions.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    with open(a.out / "damage_functions.txt", "w", encoding="utf-8") as f:
        f.write("Per PKA. Td: damage energy (TRIM's Robinson partition); TRIM-FC: TRIM full-cascade vacancies;\n"
                "repl: TRIM replacements; NRT, KP, arc-dpa, rpa: damage functions of Td (Nordlund 2018).\n"
                "Ratios are to NRT. +- is the standard error over seeds.\n\n")
        f.write(f"{'el':3s} {'Eb':>4s} {'E keV':>7s} {'Td eV':>9s} {'NRT':>8s} {'KP':>8s} {'arc':>8s} "
                f"{'TRIM-FC':>16s} {'FC/NRT':>7s} {'arc/NRT':>7s} {'repl':>8s} {'rpa':>9s}\n")
        for r in rows:
            f.write(f"{r['element']:3s} {r['Eb_eV']:4g} {r['E_keV']:7g} {r['td']:9.1f} {r['nrt']:8.2f} "
                    f"{r['kp']:8.2f} {r['arc']:8.2f} {r['trim_vac']:9.2f}+-{r['trim_vac_sem']:5.2f} "
                    f"{r['trim_vac'] / r['nrt'] if r['nrt'] else float('nan'):7.3f} "
                    f"{r['arc'] / r['nrt'] if r['nrt'] else float('nan'):7.3f} {r['trim_repl']:8.2f} "
                    f"{r['rpa']:9.2f}\n")


if __name__ == "__main__":
    main()

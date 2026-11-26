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

"""Why are ranges in crystals longer than in amorphous targets?

For several ion-crystal pairs, the mean projected range in the crystal
transport for: the amorphous target (trim2013, and the crystal transport
through an amorphous layer); the three main axes <001>, <011>, <111> at
normal incidence; a planar direction (7 deg tilt in a {100} plane); the
implanters' "random" direction (7 deg tilt, 22 deg twist); 16 random crystal
orientations; and the same directions with a static lattice and with the
electronic stopping all on the path. A scan of the tilt (twist 22 deg) shows
how the range falls off with the angle to the axis, against the critical
angle (Lindhard's psi_1, or psi_2 at low energies; summary.json's
"crystals"). A last part compares 0.5 keV B in crystalline and in
amorphous Si with the SIMS profiles of Duffy et al. (quoted by Aboy et al.,
Mater. Sci. Eng. B 124-125, 205 (2005), Fig. 1).

    python tools/crystal/range_study.py --exe build/Release/amcset.exe \
        --out build/range_study [--scale 1.0] [--only pairs,scan,boron]

Writes OUT/results.txt and OUT/runs.csv. Runs are kept in OUT/runs and not
repeated. Only the standard library is used.
"""

import argparse
import csv
import json
import math
import random
import subprocess
import sys
from pathlib import Path

# Lattice constants and Debye temperatures: Gemmell, Rev. Mod. Phys. 46, 129
# (1974), Table IV (the values AMCSET has built in).
PAIRS = [
    # name, ion, energy, structure, element, a (A), Debye (K), ions
    ("B 30 keV -> Si", "B", "30 keV", "diamond", "Si", 5.431, 543, 3000),
    ("As 100 keV -> Si", "As", "100 keV", "diamond", "Si", 5.431, 543, 1500),
    ("He 100 keV -> Si", "He", "100 keV", "diamond", "Si", 5.431, 543, 1500),
    ("Fe 20 keV -> Fe", "Fe", "20 keV", "bcc", "Fe", 2.867, 420, 2000),
    ("Cu 30 keV -> Cu", "Cu", "30 keV", "fcc", "Cu", 3.615, 315, 1500),
    ("W 30 keV -> W", "W", "30 keV", "bcc", "W", 3.165, 310, 1000),
    ("Xe 40 keV -> W", "Xe", "40 keV", "bcc", "W", 3.165, 310, 1000),
]


class Runner:
    def __init__(self, exe, out, scale):
        self.exe, self.out, self.scale = exe, out, scale
        self.lines = []
        self.rows = []

    def say(self, text=""):
        print(text, flush=True)
        self.lines.append(text)

    def run(self, name, scenario, overrides=()):
        name = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in name)
        d = self.out / "runs" / name
        summary = d / "out" / "summary.json"
        if not summary.exists():
            d.mkdir(parents=True, exist_ok=True)
            scenario = json.loads(json.dumps(scenario))
            scenario["ions"] = max(50, int(scenario["ions"] * self.scale))
            scenario.setdefault("output", {})["directory"] = str((d / "out").resolve())
            scenario["output"]["files"] = {"formats": ["json"]}
            scenario["output"]["range_3d"] = False
            scenario["output"]["exit_lists"] = False
            (d / "scenario.json").write_text(json.dumps(scenario), encoding="utf-8")
            cmd = [str(self.exe), "run", str(d / "scenario.json"), "-q"]
            for o in overrides:
                cmd += ["--set", o]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0 or not summary.exists():
                sys.exit(f"{name}: amcset failed\n{r.stdout}\n{r.stderr[-2000:]}")
        return json.loads(summary.read_text(encoding="utf-8"))


def crystal(structure, element, a, debye, surface=(0, 0, 1), reference=None, twist=0.0, euler=None):
    c = {"structure": structure, "a": f"{a} A", "sublattices": [element], "temperature": "300 K",
         "vibrations": [{"element": element, "debye_temperature": f"{debye} K"}]}
    if euler is not None:
        c["orientation"] = {"euler_deg": list(euler), "twist": f"{twist} deg"}
    else:
        o = {"surface": list(surface), "twist": f"{twist} deg"}
        if reference is not None:
            o["reference"] = list(reference)
        c["orientation"] = o
    return c


def scenario(ion, energy, ions, layer, tilt=0.0, transport="crystal", thickness="20 um"):
    layer = dict(layer)
    layer["thickness"] = thickness
    return {"ion": {"element": ion, "energy": energy, "angle": f"{tilt} deg"},
            "layers": [layer], "ions": ions, "damage": "kp",
            "physics": {"transport": transport}, "output": {"bins": 200}}


def nearest_axis(s):
    """The axis nearest the beam and its critical angle (Lindhard's psi_1,
    or psi_2 at low energies), from summary.json's crystals block."""
    try:
        a = s["crystals"][0]["axes_near_beam"][0]
        return a["uvw"], a["angle_deg"], a.get("critical_deg", a["psi1_deg"])
    except (KeyError, IndexError):
        return "", float("nan"), float("nan")


def pair_study(r, pair):
    name, ion, energy, structure, element, a, debye, ions = pair
    key = "".join(ch if ch.isalnum() else "_" for ch in name)

    def go(label, scen, ov=()):
        s = r.run(f"{key}__{label}", scen, ov)
        uvw, angle, psi1 = nearest_axis(s)
        rng = s["range"]["mean_A"]
        r.rows.append({"pair": name, "case": label, "mean_range_A": rng,
                       "straggling_A": s["range"]["straggling_A"], "vacancies_per_ion": s["per_ion"]["vacancies"],
                       "nearest_axis": uvw, "angle_deg": angle, "critical_deg": psi1})
        return s, rng, uvw, angle, psi1

    # The amorphous references: TRIM's transport, and the crystal transport
    # through an amorphous layer of the crystal's composition and density
    # (cell atoms x standard atomic mass / cell volume, as AMCSET computes it).
    cell_atoms = {"diamond": 8, "bcc": 2, "fcc": 4}[structure]
    masses = {"Si": 28.0855, "Fe": 55.845, "Cu": 63.546, "W": 183.84}
    rho = cell_atoms * masses[element] * 1.66053906660e-24 / (a * 1e-8) ** 3
    amorphous_layer = {"density": f"{rho:.6f} g/cm3", "elements": [{"element": element}]}
    _, r_trim, _, _, _ = go("amorphous_trim2013", scenario(ion, energy, ions, amorphous_layer, transport="trim2013"))
    _, r_am, _, _, _ = go("amorphous_crystal_transport", scenario(ion, energy, ions, amorphous_layer))
    ref = r_trim
    r.say(f"{name}  (amorphous: trim2013 {r_trim:.1f} A, crystal transport {r_am:.1f} A = "
          f"{r_am / r_trim:.3f})")

    def show(label, scen, ov=()):
        s, rng, uvw, angle, psi1 = go(label, scen, ov)
        r.say(f"  {label:44s} {rng:9.1f} A  x{rng / ref:5.2f}   nearest axis {uvw:9s} {angle:5.2f} deg "
              f"(critical {psi1:5.2f} deg)")
        return rng

    for label, surf in [("<001> at 0 deg", (0, 0, 1)), ("<011> at 0 deg", (0, 1, 1)), ("<111> at 0 deg", (1, 1, 1))]:
        show(label, scenario(ion, energy, ions, {"crystal": crystal(structure, element, a, debye, surf)}))
    show("(001), 7 deg tilt in a {100} plane", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye, (0, 0, 1), (1, 0, 0))}, tilt=7.0))
    show("(001), 7 deg tilt, 22 deg twist", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye, (0, 0, 1), (1, 0, 0), twist=22.0)}, tilt=7.0))
    show("(001), 7/22, static lattice", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye, (0, 0, 1), (1, 0, 0), twist=22.0)}, tilt=7.0),
        ["physics.crystal.thermal_vibrations=false"])
    show("(001), 7/22, stopping all on the path", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye, (0, 0, 1), (1, 0, 0), twist=22.0)}, tilt=7.0),
        ["physics.crystal.local_stopping_fraction=0"])
    show("<001>, static lattice", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye)}), ["physics.crystal.thermal_vibrations=false"])
    show("<001>, stopping all on the path", scenario(ion, energy, ions, {
        "crystal": crystal(structure, element, a, debye)}), ["physics.crystal.local_stopping_fraction=0"])

    # Random orientations (uniform), the same 16 for every pair.
    rng = random.Random(11)
    ranges = []
    for k in range(16):
        euler = (360.0 * rng.random(), math.degrees(math.acos(1.0 - 2.0 * rng.random())), 360.0 * rng.random())
        s, rr, uvw, angle, psi1 = go(f"random_{k:02d}", scenario(ion, energy, ions, {
            "crystal": crystal(structure, element, a, debye, euler=euler)}))
        ranges.append((rr, uvw, angle, psi1))
    ranges.sort()
    mean = sum(x[0] for x in ranges) / len(ranges)
    within = sum(1 for x in ranges if x[2] <= x[3])
    r.say(f"  16 random orientations: mean {mean:.1f} A (x{mean / ref:.2f}), least {ranges[0][0]:.1f} A "
          f"(x{ranges[0][0] / ref:.2f}), median {0.5 * (ranges[7][0] + ranges[8][0]):.1f} A, most "
          f"{ranges[-1][0]:.1f} A (x{ranges[-1][0] / ref:.2f}); the beam within the critical angle of an axis "
          f"(indices <= 2) in {within} of 16")
    r.say()


def tilt_scan(r, pair, tilts):
    name, ion, energy, structure, element, a, debye, ions = pair
    key = "".join(ch if ch.isalnum() else "_" for ch in name)
    r.say(f"{name}: tilt from <001> (twist 22 deg, so that no low-index plane contains the beam)")
    for tilt in tilts:
        s = r.run(f"scan_{key}__{tilt}", scenario(ion, energy, ions, {
            "crystal": crystal(structure, element, a, debye, (0, 0, 1), (1, 0, 0), twist=22.0)}, tilt=tilt))
        uvw, angle, psi1 = nearest_axis(s)
        r.say(f"  tilt {tilt:5.1f} deg  range {s['range']['mean_A']:8.1f} A   nearest axis {uvw} at {angle:.2f} deg "
              f"(critical {psi1:.2f})")
        r.rows.append({"pair": name, "case": f"scan tilt {tilt}", "mean_range_A": s["range"]["mean_A"],
                       "straggling_A": s["range"]["straggling_A"], "vacancies_per_ion": s["per_ion"]["vacancies"],
                       "nearest_axis": uvw, "angle_deg": angle, "critical_deg": psi1})
    r.say()


def profile(results_path):
    t = None
    for tb in json.loads(results_path.read_text(encoding="utf-8"))["tables"]:
        if tb["key"] == "range":
            t = tb["data"]
    width = t["depth"][1] - t["depth"][0]
    return [d - 0.5 * width for d in t["depth"]], t["ions"]


def depth_at(depth, p, level):
    s = [sum(p[max(0, i - 1):i + 2]) / len(p[max(0, i - 1):i + 2]) for i in range(len(p))]
    peak = max(range(len(s)), key=lambda i: s[i])
    last = float("nan")
    for i in range(peak, len(s) - 1):
        if s[i] >= level > s[i + 1]:
            last = depth[i]
    return last


def boron(r):
    # Duffy et al. (Aboy et al. 2005, Fig. 1, as implanted; 0.5 keV, 1e15
    # cm^-2; read off the figure, about 0.3 nm): depth (nm) at 1e20, 1e19 and
    # 3e18 cm^-3 in crystalline and in preamorphized Si (amorphous to 29 nm).
    measured = {"crystalline": (8.2, 15.0, float("nan")), "preamorphized": (8.3, 11.3, 14.5)}
    dose = 1e15
    r.say("0.5 keV B, 1e15 cm-2, into crystalline and preamorphized Si: depth (nm) where the concentration falls "
          "through")
    r.say("  (measured: Duffy et al. as quoted by Aboy et al., MSE B 124-125, 205 (2005), Fig. 1; tilt not stated)")
    r.say(f"  {'':46s} {'1e20':>6s} {'1e19':>6s} {'3e18':>6s} cm-3")
    for k, v in measured.items():
        r.say(f"  {'measured, ' + k:46s} " + " ".join(f"{x:6.1f}" for x in v))
    base = {"ion": {"element": "B", "energy": "0.5 keV", "angle": "0 deg"}, "ions": 400000, "damage": "kp",
            "output": {"depth_window": ["0 A", "400 A"], "bins": 160}}
    si = {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"], "temperature": "300 K",
          "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0]}}
    cases = [("amorphous (trim2013)", "trim2013", None, 0.0, []),
             ("amorphous (crystal transport)", "crystal", None, 0.0, []),
             ("crystal <001>, 0 deg", "crystal", si, 0.0, []),
             ("crystal, 7 deg tilt, 22 deg twist", "crystal", si, 7.0, ["layers.0.crystal.orientation.twist=22 deg"]),
             ("crystal <001>, 0 deg, stopping all on the path", "crystal", si, 0.0,
              ["physics.crystal.local_stopping_fraction=0"])]
    for label, transport, xt, tilt, ov in cases:
        scen = json.loads(json.dumps(base))
        scen["ion"]["angle"] = f"{tilt} deg"
        scen["physics"] = {"transport": transport}
        if xt is None:
            scen["layers"] = [{"thickness": "3000 A", "density": "2.3291 g/cm3", "elements": [{"element": "Si"}]}]
        else:
            scen["layers"] = [{"thickness": "3000 A", "crystal": xt}]
        name = "boron05_" + "".join(ch if ch.isalnum() else "_" for ch in label)
        r.run(name, scen, ov)
        d, p = profile(r.out / "runs" / name / "out" / "results.json")
        c = [v * dose for v in p]
        r.say(f"  {label:46s} " + " ".join(f"{depth_at(d, c, lv) / 10:6.1f}" for lv in (1e20, 1e19, 3e18)))
    r.say()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--only", default="pairs,scan,boron")
    a = ap.parse_args()
    r = Runner(a.exe, a.out, a.scale)
    version = subprocess.run([str(a.exe), "version"], capture_output=True, text=True).stdout.strip()
    r.say(f"Ranges in crystals and amorphous targets, {version}, ions x {a.scale:g}")
    r.say("SRIM-2013 stopping; crystals at 300 K with Gemmell's Debye temperatures; quick (Kinchin-Pease) damage;")
    r.say("x = the mean range over that of the amorphous target in TRIM's transport (trim2013).")
    r.say()
    parts = a.only.split(",")
    if "pairs" in parts:
        for p in PAIRS:
            pair_study(r, p)
    if "scan" in parts:
        tilt_scan(r, PAIRS[0], [0, 0.5, 1, 1.5, 2, 3, 4, 5, 7, 10, 15])
        tilt_scan(r, PAIRS[3], [0, 1, 2, 3, 4, 5, 7, 10, 15, 20])
    if "boron" in parts:
        boron(r)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "results.txt").write_text("\n".join(r.lines) + "\n", encoding="utf-8")
    if r.rows:
        with (a.out / "runs.csv").open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(r.rows[0].keys()))
            w.writeheader()
            w.writerows(r.rows)


if __name__ == "__main__":
    main()

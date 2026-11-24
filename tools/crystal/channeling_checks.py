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

"""Checks of the crystal transport (physics.transport: crystal).

Runs the cases whose numbers docs/CRYSTAL_DESIGN.md and docs/CLI.md quote and
writes them to OUT/results.txt:

  directions   30 keV B into Si: <001>, <110>, 7 deg tilt / 22 deg twist, a
               100 A oxide on top, against the amorphous (trim2013) target
  parameters   the same along <001>: temperature, static lattice, how the
               electronic stopping is divided, simultaneity distance,
               maximum impact parameter
  dip          1 MeV He through 600 A of W along <001>: knock-on atoms per
               ion against the tilt angle (Lindhard's psi_1 = 2.10 deg)
  straight     1 MeV H through 3000 A of Si in a random direction: the
               electronic loss with the stopping on the path, in the
               collisions, and in TRIM (the normalisation)
  random       averages over random crystal orientations against the
               amorphous target (He in W, Fe in Fe, B in Si)
  cutoff       full cascades with TRIM's 1 eV cut-off and with the cut-off at
               Ed, against NRT

    python tools/crystal/channeling_checks.py --exe build/Release/amcset.exe \
        --out build/crystal_checks [--only dip,cutoff] [--scale 1.0]

Only the standard library is used. --scale multiplies the numbers of ions.
These are checks of trends and of internal consistency, not comparisons with
measurements: the transport's electronic stopping in channels is not
calibrated (docs/CRYSTAL_DESIGN.md).
"""

import argparse
import json
import math
import random
import subprocess
import sys
from pathlib import Path

SI = {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
      "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0]},
      "vibrations": [{"element": "Si", "debye_temperature": "645 K"}]}
W = {"structure": "bcc", "a": "3.165 A", "sublattices": ["W"],
     "orientation": {"surface": [0, 0, 1], "reference": [1, 0, 0]},
     "vibrations": [{"element": "W", "debye_temperature": "310 K"}]}
FE = {"structure": "bcc", "a": "2.8665 A", "sublattices": ["Fe"],
      "vibrations": [{"element": "Fe", "debye_temperature": "470 K"}]}


class Runner:
    def __init__(self, exe, out, scale):
        self.exe, self.out, self.scale = exe, out, scale
        self.lines = []

    def say(self, text=""):
        print(text, flush=True)
        self.lines.append(text)

    def run(self, name, scenario, overrides=()):
        d = self.out / "runs" / name
        d.mkdir(parents=True, exist_ok=True)
        scenario = json.loads(json.dumps(scenario))
        scenario["ions"] = max(20, int(scenario["ions"] * self.scale))
        scenario.setdefault("output", {})["directory"] = str((d / "out").resolve())
        scenario["output"]["files"] = {"formats": ["json"]}
        (d / "scenario.json").write_text(json.dumps(scenario), encoding="utf-8")
        cmd = [str(self.exe), "run", str(d / "scenario.json"), "-q"]
        for o in overrides:
            cmd += ["--set", o]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"{name}: amcset failed\n{r.stdout}\n{r.stderr[-2000:]}")
        return json.loads((d / "out" / "summary.json").read_text(encoding="utf-8"))


def b_si(ions=2000, damage="kp"):
    return {"ion": {"element": "B", "energy": "30 keV"}, "layers": [{"thickness": "4 um", "crystal": SI}],
            "ions": ions, "damage": damage, "physics": {"transport": "crystal"}}


TILTED = ["ion.angle=7 deg", "layers.0.crystal.orientation.twist=22 deg"]


def row(r, label, s):
    r.say(f"  {label:34s} range {s['range']['mean_A']:8.1f} A  straggling {s['range']['straggling_A']:7.1f}  "
          f"lateral {s['range']['lateral_A']:6.1f}  vacancies/ion {s['per_ion']['vacancies']:7.2f}")


def directions(r):
    r.say("30 keV B into Si (300 K, Debye temperature 645 K), quick damage, 2000 ions")
    row(r, "<001>", r.run("dir_001", b_si()))
    row(r, "<110>", r.run("dir_110", b_si(), ["layers.0.crystal.orientation.surface=[1,1,0]",
                                              "layers.0.crystal.orientation.reference=[0,0,1]"]))
    row(r, "7 deg tilt, 22 deg twist", r.run("dir_tilt", b_si(), TILTED))
    oxide = b_si()
    oxide["layers"].insert(0, {"name": "SiO2", "thickness": "100 A", "density": "2.2 g/cm3",
                               "elements": [{"element": "Si", "stoichiometry": 1},
                                            {"element": "O", "stoichiometry": 2}]})
    row(r, "<001> under 100 A of SiO2", r.run("dir_oxide", oxide))
    row(r, "amorphous (crystal transport)", r.run("dir_amorphous_ct", {
        "ion": {"element": "B", "energy": "30 keV"}, "layers": [{"thickness": "4 um", "elements": [{"element": "Si"}]}],
        "ions": 2000, "damage": "kp", "physics": {"transport": "crystal"}}))
    row(r, "amorphous (trim2013)", r.run("dir_trim", b_si(), ["physics.transport=trim2013"]))
    r.say()


def parameters(r):
    r.say("30 keV B into Si <001>, 2000 ions: what the channeled range depends on")
    for label, ov in [("20 K", ["layers.0.crystal.temperature=20 K"]),
                      ("300 K", []),
                      ("1200 K", ["layers.0.crystal.temperature=1200 K"]),
                      ("static lattice", ["physics.crystal.thermal_vibrations=false"]),
                      ("stopping all on the path", ["physics.crystal.local_stopping_fraction=0"]),
                      ("stopping all in the collisions", ["physics.crystal.local_stopping_fraction=1"]),
                      ("sequential collisions (0 A)", ["physics.crystal.simultaneous_distance=0 A"]),
                      ("simultaneity distance 1.5 A", ["physics.crystal.simultaneous_distance=1.5 A"]),
                      ("maximum impact parameter 3 A", ["physics.crystal.max_impact_parameter=3 A"])]:
        row(r, label, r.run("par_" + "".join(c if c.isalnum() else "_" for c in label), b_si(), ov))
    r.say("The same 7 deg off, 22 deg twist (a random direction):")
    for label, ov in [("stopping all on the path", ["physics.crystal.local_stopping_fraction=0"]),
                      ("stopping all in the collisions", ["physics.crystal.local_stopping_fraction=1"]),
                      ("maximum impact parameter 3 A", ["physics.crystal.max_impact_parameter=3 A"])]:
        row(r, label, r.run("part_" + "".join(c if c.isalnum() else "_" for c in label), b_si(), TILTED + ov))
    r.say()


def dip(r):
    psi1 = math.degrees(math.sqrt(2 * 2 * 74 * 14.4 / (1e6 * 3.165)))
    r.say(f"1 MeV He through 600 A of W along <001> (300 K): atoms displaced by the ion, per ion; "
          f"psi_1 = {psi1:.2f} deg; 6000 ions per angle; twist 17 deg")
    base = {"ion": {"element": "He", "energy": "1 MeV"}, "layers": [{"thickness": "600 A", "crystal": W}],
            "ions": 6000, "damage": "kp", "physics": {"transport": "crystal"}}
    ref = r.run("dip_trim", base, ["physics.transport=trim2013"])["per_ion"]["knock_ons"]
    r.say(f"  amorphous (trim2013): {ref:.4f}")
    for tilt in [0.0, 0.5, 1.0, 1.25, 1.5, 1.75, 2.1, 2.5, 3.0, 4.0, 6.0]:
        s = r.run(f"dip_{tilt}", base, [f"ion.angle={tilt} deg", "layers.0.crystal.orientation.twist=17 deg"])
        k = s["per_ion"]["knock_ons"]
        r.say(f"  tilt {tilt:5.2f} deg  {k:8.4f}  = {k / ref:6.3f} of the amorphous value")
    r.say()


def straight(r):
    r.say("1 MeV H through 3000 A of Si, 11 deg tilt / 27 deg twist, no straggling: energy lost to electrons (%)")
    base = {"ion": {"element": "H", "energy": "1 MeV", "angle": "11 deg"},
            "layers": [{"thickness": "3000 A", "crystal": SI}], "ions": 400, "damage": "kp",
            "physics": {"transport": "crystal", "electronic_straggling": False}}
    twist = ["layers.0.crystal.orientation.twist=27 deg"]
    for label, ov in [("all on the path", ["physics.crystal.local_stopping_fraction=0"]),
                      ("half and half (default)", []),
                      ("all in the collisions", ["physics.crystal.local_stopping_fraction=1"]),
                      ("all in the collisions, s_max 3.5 A", ["physics.crystal.local_stopping_fraction=1",
                                                              "physics.crystal.max_impact_parameter=3.5 A"]),
                      ("trim2013 (amorphous)", ["physics.transport=trim2013"])]:
        s = r.run("str_" + "".join(c if c.isalnum() else "_" for c in label), base, twist + ov)
        r.say(f"  {label:36s} {s['energy_percent']['ionization_ion']:.4f}")
    r.say()


def random_orientations(r):
    r.say("Averages over 16 random crystal orientations (quick damage) against the amorphous target")
    rng = random.Random(20261124)
    for label, ion, energy, layer, cell, ions, key in [
            ("1 MeV He, 600 A W", "He", "1 MeV", {"thickness": "600 A"}, W, 800, "knock_ons"),
            ("20 keV Fe in Fe", "Fe", "20 keV", {"thickness": "3000 A", "elements": [{"element": "Fe", "Ed": "40 eV", "Eb": "0 eV"}]}, FE, 150, "vacancies"),
            ("30 keV B in Si", "B", "30 keV", {"thickness": "2 um"}, SI, 300, "vacancies"),
            ("200 keV Fe in Fe", "Fe", "200 keV", {"thickness": "1 um", "elements": [{"element": "Fe", "Ed": "40 eV", "Eb": "0 eV"}]}, FE, 100, "vacancies")]:
        values, ranges = [], []
        tag = "".join(c if c.isalnum() else "_" for c in label)
        for k in range(16):
            cr = dict(cell)
            cr["orientation"] = {"euler_deg": [360 * rng.random(), math.degrees(math.acos(1 - 2 * rng.random())),
                                               360 * rng.random()]}
            lay = dict(layer)
            lay["crystal"] = cr
            s = r.run(f"rnd_{tag}_{k}", {"ion": {"element": ion, "energy": energy}, "layers": [lay], "ions": ions,
                                         "seed": 100 + k, "damage": "kp", "physics": {"transport": "crystal"}})
            values.append(s["per_ion"][key])
            ranges.append(s["range"]["mean_A"])
        lay = dict(layer)
        lay["crystal"] = dict(cell)
        t = r.run(f"rnd_{tag}_trim", {"ion": {"element": ion, "energy": energy}, "layers": [lay], "ions": ions * 16,
                                      "damage": "kp", "physics": {"transport": "trim2013"}})

        def stats(v):
            m = sum(v) / len(v)
            sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
            return m, sd / math.sqrt(len(v))
        m, se = stats(values)
        r.say(f"  {label:20s} {key}/ion {m:8.4g} +- {se:.2g} (amorphous {t['per_ion'][key]:.4g})")
        if max(ranges) > 0:
            mr, ser = stats(ranges)
            r.say(f"  {'':20s} range {mr:7.1f} +- {ser:5.1f} A, smallest {min(ranges):.0f}, largest {max(ranges):.0f} "
                  f"(amorphous {t['range']['mean_A']:.1f})")
    r.say()


def cutoff(r):
    r.say("30 keV B into Si (Ed = 15 eV, Eb = 0), 7 deg tilt / 22 deg twist, 2000 ions: vacancies per ion")
    base = b_si(damage="full")
    base["layers"][0]["elements"] = [{"element": "Si", "Ed": "15 eV", "Eb": "0 eV"}]
    s = r.run("cut_1eV", base, TILTED)
    r.say(f"  crystal, full cascades, 1 eV cut-off (TRIM's)   {s['per_ion']['vacancies']:7.1f}  replacements {s['per_ion']['replacements']:6.1f}")
    s = r.run("cut_Ed", base, TILTED + ["physics.recoil_stop_energy=15 eV", "physics.ion_stop_energy=15 eV"])
    r.say(f"  crystal, full cascades, cut-off at Ed            {s['per_ion']['vacancies']:7.1f}  replacements {s['per_ion']['replacements']:6.1f}")
    s = r.run("cut_trim_full", base, ["physics.transport=trim2013"])
    r.say(f"  trim2013 (amorphous), full cascades              {s['per_ion']['vacancies']:7.1f}  replacements {s['per_ion']['replacements']:6.1f}")
    s = r.run("cut_trim_kp", base, ["physics.transport=trim2013", "damage=kp"])
    r.say(f"  trim2013 quick calculation (NRT)                 {s['per_ion']['vacancies']:7.1f}")
    r.say()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only", default="directions,parameters,dip,straight,random,cutoff")
    ap.add_argument("--scale", type=float, default=1.0)
    a = ap.parse_args()
    r = Runner(a.exe, a.out, a.scale)
    version = subprocess.run([str(a.exe), "version"], capture_output=True, text=True).stdout.strip()
    r.say(f"Crystal transport checks, {version}, ions x {a.scale:g}")
    r.say("Stopping: the SRIM-2013 tables that come with AMCSET. Trends and consistency, not measurements.")
    r.say()
    parts = {"directions": directions, "parameters": parameters, "dip": dip, "straight": straight,
             "random": random_orientations, "cutoff": cutoff}
    for name in a.only.split(","):
        if name not in parts:
            sys.exit(f"unknown part '{name}' (one of {', '.join(parts)})")
        parts[name](r)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "results.txt").write_text("\n".join(r.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

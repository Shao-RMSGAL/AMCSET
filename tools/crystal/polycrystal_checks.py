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

"""Checks of polycrystals in the crystal transport: trends, not measurements.

30 keV B into silicon along the normal of (0 0 1), a channel of the single
crystal: the mean range and the vacancies per ion for the single crystal,
the amorphous target, and grains of several sizes, textures and boundary
widths; and the same grains read back from a grain file.

    python tools/crystal/polycrystal_checks.py --exe build/Release/amcset.exe \
        --out build/polycrystal_checks [--scale 1.0]

Writes OUT/results.txt (the source of the numbers in docs/CLI.md "Crystals").
Only the standard library is used. Runs are kept in OUT/runs and not
repeated.
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
      "temperature": "300 K", "vibrations": [{"element": "Si", "debye_temperature": "543 K"}]}


class Runner:
    def __init__(self, exe, out, scale):
        self.exe, self.out, self.scale = exe, out, scale
        self.lines = []

    def say(self, text=""):
        print(text, flush=True)
        self.lines.append(text)

    def run(self, name, scenario, overrides=()):
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


def scenario(polycrystal=None, ions=4000, spot="1 um"):
    crystal = json.loads(json.dumps(SI))
    if polycrystal is not None:
        crystal["polycrystal"] = polycrystal
    return {"ion": {"element": "B", "energy": "30 keV", "angle": "0 deg",
                    "spot": {"shape": "rectangle", "width": spot, "height": spot}},
            "layers": [{"thickness": "3 um", "crystal": crystal}],
            "ions": ions, "damage": "kp", "physics": {"transport": "crystal"}}


def row(r, label, s):
    r.say(f"  {label:46s} range {s['range']['mean_A']:8.1f} A   straggling {s['range']['straggling_A']:7.1f} A   "
          f"vacancies/ion {s['per_ion']['vacancies']:6.1f}   {s['elapsed_seconds']:6.2f} s")
    return s["range"]["mean_A"]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--scale", type=float, default=1.0)
    a = ap.parse_args()
    r = Runner(a.exe, a.out, a.scale)
    version = subprocess.run([str(a.exe), "version"], capture_output=True, text=True).stdout.strip()
    r.say(f"Polycrystals in the crystal transport, {version}, ions x {a.scale:g}")
    r.say("30 keV B into Si along the normal of (0 0 1), 300 K; 4000 ions each; SRIM-2013 stopping.")
    r.say()

    r.say("The limits:")
    row(r, "single crystal, <001>", r.run("single", scenario()))
    row(r, "single crystal, 7 deg tilt, 22 deg twist", r.run(
        "single_tilted", scenario(), ["ion.angle=7 deg", "layers.0.crystal.orientation.twist=22 deg"]))
    amorphous = scenario()
    del amorphous["layers"][0]["crystal"]
    amorphous["layers"][0].update({"density": "2.3291 g/cm3", "elements": [{"element": "Si"}]})
    row(r, "amorphous (crystal transport, no crystal)", r.run("amorphous", amorphous))
    r.say()

    r.say("Grains of random orientations (the beam over one period of the grains):")
    for size in (50, 100, 200, 500, 1000, 3000):
        p = {"source": "voronoi", "grain_size": f"{size} A", "texture": "random"}
        row(r, f"grain size {size} A", r.run(f"random_{size}", scenario(p, spot=f"{8 * size} A")))
    p = {"source": "voronoi", "grain_size": "200 A", "texture": "random", "columnar": True}
    row(r, "columns of 200 A", r.run("random_columns_200", scenario(p, spot="1600 A")))
    p = {"source": "voronoi", "grain_size": "200 A", "texture": "random", "seed": 2}
    row(r, "grain size 200 A, other grains (seed 2)", r.run("random_200_seed2", scenario(p, spot="1600 A")))
    r.say()

    r.say("Textures, grains of 500 A:")
    for label, extra in [("fiber, <001> along the normal", {"texture": "fiber"}),
                         ("fiber, 2 deg spread", {"texture": "fiber", "spread": "2 deg"}),
                         ("fiber, 5 deg spread", {"texture": "fiber", "spread": "5 deg"}),
                         ("mosaic, no spread (lattices not in registry)", {"texture": "mosaic"}),
                         ("mosaic, 0.5 deg spread", {"texture": "mosaic", "spread": "0.5 deg"}),
                         ("mosaic, 2 deg spread", {"texture": "mosaic", "spread": "2 deg"})]:
        p = {"source": "voronoi", "grain_size": "500 A"}
        p.update(extra)
        name = "tex_" + "".join(c if c.isalnum() else "_" for c in label)
        row(r, label, r.run(name, scenario(p, spot="4000 A")))
    r.say()

    r.say("Disordered grain boundaries (treated as amorphous), fiber texture, grains of 100 A:")
    for width in (0, 5, 10, 20, 40):
        p = {"source": "voronoi", "grain_size": "100 A", "texture": "fiber", "boundary_width": f"{width} A"}
        row(r, f"boundary width {width} A", r.run(f"width_{width}", scenario(p, spot="800 A")))
    r.say()

    # The same kind of grains from a grain file: seeds and Bunge angles
    # written here, read by AMCSET.
    r.say("Grains from a grain file (400 seeds in 1600 x 1600 x 3400 A, random orientations):")
    rng = random.Random(5)
    lines = ["# written by tools/crystal/polycrystal_checks.py", "unit A", "period 1600 1600"]
    for _ in range(400):
        lines.append("grain %.3f %.3f %.3f %.4f %.4f %.4f" % (
            -200 + 3400 * rng.random(), 1600 * rng.random(), 1600 * rng.random(),
            360 * rng.random(), math.degrees(math.acos(1 - 2 * rng.random())), 360 * rng.random()))
    a.out.mkdir(parents=True, exist_ok=True)
    grain_file = a.out / "grains.txt"
    grain_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
    p = {"file": str(grain_file.resolve())}
    file_scenario = scenario(p, spot="1600 A")
    file_scenario["layers"][0]["thickness"] = "3000 A"
    row(r, "grain file, 3000 A of Si", r.run("grain_file", file_scenario))
    generated = scenario({"source": "voronoi", "grain_size": "254 A", "texture": "random",
                          "period": ["1600 A", "1600 A"]}, spot="1600 A")
    generated["layers"][0]["thickness"] = "3000 A"
    row(r, "generated, the same number of grains", r.run("grain_file_like", generated))

    (a.out / "results.txt").write_text("\n".join(r.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

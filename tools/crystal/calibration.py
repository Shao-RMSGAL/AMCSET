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

"""The crystal transport against published measurements.

  dip     1 MeV He along <110> of Ge at 20 C: the channeling dip (atoms
          displaced by the ion behind the surface, against the tilt angle).
          Measured half-angle 0.95 deg and minimum yield 0.03 (Picraux et al.,
          as quoted by Gemmell, Rev. Mod. Phys. 46, 129 (1974), Appendix A,
          with Barrett's formulas: 0.93 deg, 0.025).
  boron   Boron implanted into (100) Si, as-implanted SIMS profiles, for
          several divisions of the electronic stopping between the path and
          the collisions (physics.crystal.local_stopping_fraction):
            - 15 keV, 3e15 cm-2, 7 deg tilt: the moments of the profile
              (Jung, Trans. Electr. Electron. Mater. 11, 120 (2010), Table 2);
            - 30 keV, 1.5e14 cm-2: peak and tail (Shao et al., Mater. Sci.
              Eng. R 42, 65 (2003), Fig. 33; tilt not stated);
            - 0.5 keV, 1e14 and 1e15 cm-2: the tail (Schroer et al., Mater.
              Sci. Eng. B 71, 219 (2000), Fig. 1; tilt not stated).
          The profile values were read off the published figures (a few per
          cent of the depth); the moments are tabulated in the paper.

    python tools/crystal/calibration.py --exe build/Release/amcset.exe \
        --out build/crystal_calibration [--only dip,boron] [--scale 1.0] [--fractions 0,0.25,0.5,0.75,1]

Writes OUT/results.txt. Only the standard library is used. The simulation is
of a perfect crystal (every ion sees an undamaged lattice): the measurements
at 1e15 cm-2 and above include the dechanneling by the damage of the earlier
ions, which it does not.
"""

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path


class Runner:
    def __init__(self, exe, out, scale):
        self.exe, self.out, self.scale = exe, out, scale
        self.lines = []

    def say(self, text=""):
        print(text, flush=True)
        self.lines.append(text)

    def run(self, name, scenario, overrides=()):
        d = self.out / "runs" / name
        results = d / "out" / "results.json"
        if not results.exists():
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
            if r.returncode != 0 or not results.exists():
                sys.exit(f"{name}: amcset failed\n{r.stdout}\n{r.stderr[-2000:]}")
        return json.loads(results.read_text(encoding="utf-8")), json.loads((d / "out" / "summary.json").read_text(encoding="utf-8"))


def table(results, key):
    for t in results["tables"]:
        if t["key"] == key:
            return t["data"]  # column key -> values
    raise KeyError(key)


# ---------------------------------------------------------------------------------------
# The channeling dip
# ---------------------------------------------------------------------------------------

def dip(r):
    a = 5.657  # Ge lattice constant (Gemmell, Table IV), Angstrom
    psi1 = math.degrees(math.sqrt(2 * 2 * 32 * 14.4 / (1e6 * a / math.sqrt(2))))
    r.say(f"1 MeV He along <110> of Ge at 293 K (Debye temperature 290 K): atoms displaced by the ion at depths of "
          f"200-1200 A, relative to the amorphous target; psi_1 = {psi1:.2f} deg")
    r.say("Measured (Picraux et al., quoted by Gemmell 1974, App. A): half-angle 0.95 deg, minimum yield 0.03; "
          "Barrett's formulas: 0.93 deg, 0.025")
    crystal = {"structure": "diamond", "a": f"{a} A", "sublattices": ["Ge"],
               "orientation": {"surface": [1, 1, 0], "reference": [0, 0, 1], "twist": "17 deg"},
               "temperature": "293 K", "vibrations": [{"element": "Ge", "debye_temperature": "290 K"}]}
    base = {"ion": {"element": "He", "energy": "1 MeV"}, "layers": [{"thickness": "1500 A", "crystal": crystal}],
            "ions": 40000, "damage": "kp", "physics": {"transport": "crystal"},
            "output": {"bins": 150}}

    def yield_of(name, overrides):
        res, _ = r.run(name, base, overrides)
        v = table(res, "vacancies")
        width = v["depth"][1] - v["depth"][0]
        return sum(k * width for d, k in zip(v["depth"], v["by_ions"]) if 200.0 < d <= 1200.0)

    random = yield_of("dip_amorphous", ["physics.transport=trim2013"])
    r.say(f"  amorphous (trim2013): {random:.4f} per ion")
    angles = [0.0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.4, 1.7, 2.0, 2.5]
    chi = []
    for tilt in angles:
        y = yield_of(f"dip_{tilt}", [f"ion.angle={tilt} deg"])
        chi.append(y / random)
        r.say(f"  tilt {tilt:4.2f} deg   {y:8.5f}   chi = {chi[-1]:6.3f}")
    chi_min = chi[0]
    half = 0.5 * (1.0 + chi_min)
    width = float("nan")
    for i in range(1, len(angles)):
        if chi[i - 1] < half <= chi[i]:
            width = angles[i - 1] + (half - chi[i - 1]) / (chi[i] - chi[i - 1]) * (angles[i] - angles[i - 1])
            break
    r.say(f"  minimum yield {chi_min:.3f} (measured 0.03); half-angle at chi = {half:.2f}: {width:.2f} deg "
          f"(measured 0.95) = {width / psi1:.2f} psi_1")
    r.say()


# ---------------------------------------------------------------------------------------
# Boron in silicon
# ---------------------------------------------------------------------------------------

SI = {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
      "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0]},
      "temperature": "300 K", "vibrations": [{"element": "Si", "debye_temperature": "543 K"}]}


def boron_scenario(energy, tilt, twist, ions, window_A, bins, thickness):
    crystal = json.loads(json.dumps(SI))
    crystal["orientation"]["twist"] = f"{twist} deg"
    return {"ion": {"element": "B", "energy": energy, "angle": f"{tilt} deg"},
            "layers": [{"thickness": thickness, "crystal": crystal}],
            "ions": ions, "damage": "kp", "physics": {"transport": "crystal"},
            "output": {"depth_window": ["0 A", f"{window_A} A"], "bins": bins}}


def profile(results):
    """Depth (bin centre, A) and ions per cm per ion (= concentration / dose)."""
    t = table(results, "range")
    width = t["depth"][1] - t["depth"][0]
    return [d - 0.5 * width for d in t["depth"]], t["ions"], width


def moments(depth, p, width, upto_A):
    xs = [(x, v * width) for x, v in zip(depth, p) if x <= upto_A]
    n = sum(w for _, w in xs)
    mean = sum(x * w for x, w in xs) / n
    m2 = sum((x - mean) ** 2 * w for x, w in xs) / n
    m3 = sum((x - mean) ** 3 * w for x, w in xs) / n
    m4 = sum((x - mean) ** 4 * w for x, w in xs) / n
    sd = math.sqrt(m2)
    return mean, sd, m3 / sd ** 3, m4 / m2 ** 2


def depth_at(depth, p, level):
    """The deepest depth (A) at which the profile (smoothed over 3 bins) falls through \\p level."""
    s = [sum(p[max(0, i - 1):i + 2]) / len(p[max(0, i - 1):i + 2]) for i in range(len(p))]
    peak = max(range(len(s)), key=lambda i: s[i])
    last = float("nan")
    for i in range(peak, len(s) - 1):
        if s[i] >= level > s[i + 1] and s[i + 1] > 0:
            f = (math.log(s[i]) - math.log(level)) / (math.log(s[i]) - math.log(s[i + 1]))
            last = depth[i] + f * (depth[i + 1] - depth[i])
        elif s[i] >= level > s[i + 1]:
            last = depth[i]
    return last


def boron(r, fractions):
    # --- 15 keV, 7 deg tilt: the moments ------------------------------------------------
    r.say("15 keV B into (100) Si, 7 deg tilt, 3e15 cm-2 (Jung 2010, Table 2). Moments of the profile to 450 nm:")
    r.say("  measured (SIMS)              Rp  64.2 nm   dRp 31.8 nm   skewness 0.98   kurtosis 5.56")
    r.say("  UT-Marlowe (the paper)       Rp  59.7 nm   dRp 30.1 nm   skewness 0.72   kurtosis 4.62")
    for label, tilt, twist, ov in ([("trim2013 (amorphous)", 7, 22, ["physics.transport=trim2013"])] +
                                   [(f"crystal, fraction {f:g}, twist 22", 7, 22, [f"physics.crystal.local_stopping_fraction={f}"])
                                    for f in fractions] +
                                   [("crystal, fraction 0.5, twist 0", 7, 0, []),
                                    ("crystal, fraction 0.5, twist 45", 7, 45, []),
                                    ("crystal, fraction 0.5, 0 deg tilt", 0, 0, [])]):
        name = "b15_" + "".join(c if c.isalnum() else "_" for c in label)
        res, _ = r.run(name, boron_scenario("15 keV", tilt, twist, 100000, 6000, 300, "2 um"), ov)
        d, p, w = profile(res)
        m = moments(d, p, w, 4500.0)
        r.say(f"  {label:28s} Rp {m[0] / 10:5.1f} nm   dRp {m[1] / 10:4.1f} nm   skewness {m[2]:4.2f}   kurtosis {m[3]:4.2f}")
    r.say()

    # --- 30 keV: peak and tail -------------------------------------------------------------
    dose = 1.5e14
    r.say("30 keV B into Si, 1.5e14 cm-2, as implanted (Shao 2003, Fig. 33, in regrown float-zone Si; tilt not stated).")
    r.say("  read off the figure:         peak 1e19 cm-3 at about 105 nm;  1e18 at about 210 nm;  1e17 at about 290 nm")
    for label, tilt, twist, ov in ([("trim2013 (amorphous)", 7, 22, ["physics.transport=trim2013"])] +
                                   [(f"crystal 7/22, fraction {f:g}", 7, 22, [f"physics.crystal.local_stopping_fraction={f}"])
                                    for f in fractions] +
                                   [("crystal 0 deg, fraction 0.5", 0, 0, [])]):
        name = "b30_" + "".join(c if c.isalnum() else "_" for c in label)
        res, _ = r.run(name, boron_scenario("30 keV", tilt, twist, 100000, 8000, 320, "3 um"), ov)
        d, p, w = profile(res)
        c = [v * dose for v in p]  # cm-3
        peak = max(range(len(c)), key=lambda i: sum(c[max(0, i - 1):i + 2]))
        r.say(f"  {label:28s} peak {sum(c[peak - 1:peak + 2]) / 3:8.2e} at {d[peak] / 10:5.0f} nm;  1e18 at "
              f"{depth_at(d, c, 1e18) / 10:5.0f} nm;  1e17 at {depth_at(d, c, 1e17) / 10:5.0f} nm")
    r.say()

    # --- 0.5 keV: the tail -------------------------------------------------------------------
    r.say("0.5 keV B into (100) Si, oxide stripped, as implanted (Schroer 2000, Fig. 1; tilt not stated).")
    r.say("Depth (nm) at which the concentration falls through the level, as concentration / dose (1/cm):")
    r.say("  level (1/cm)                     1e5     1e4     1e3     1e2")
    r.say("  measured, 1e14 cm-2             14.2    25.2    35      57   (1e19 ... 1e16 cm-3)")
    r.say("  measured, 1e15 cm-2             10.1    15.9    29      49   (1e20 ... 1e17 cm-3)")
    for label, tilt, ov in ([("trim2013 (amorphous)", 0, ["physics.transport=trim2013"])] +
                            [(f"crystal 0 deg, fraction {f:g}", 0, [f"physics.crystal.local_stopping_fraction={f}"])
                             for f in fractions] +
                            [("crystal 7/22, fraction 0.5", 7, ["layers.0.crystal.orientation.twist=22 deg"]),
                             ("crystal 0 deg, static lattice", 0, ["physics.crystal.thermal_vibrations=false"])]):
        name = "b05_" + "".join(c if c.isalnum() else "_" for c in label)
        res, _ = r.run(name, boron_scenario("0.5 keV", tilt, 0, 2000000, 1000, 200, "3000 A"), ov)
        d, p, w = profile(res)
        r.say(f"  {label:30s} " + "  ".join(f"{depth_at(d, p, lv) / 10:6.1f}" for lv in (1e5, 1e4, 1e3, 1e2)))
    r.say()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only", default="dip,boron")
    ap.add_argument("--scale", type=float, default=1.0)
    ap.add_argument("--fractions", default="0,0.25,0.5,0.75,1")
    a = ap.parse_args()
    r = Runner(a.exe, a.out, a.scale)
    version = subprocess.run([str(a.exe), "version"], capture_output=True, text=True).stdout.strip()
    r.say(f"Crystal transport against measurements, {version}, ions x {a.scale:g}")
    r.say("Stopping: the SRIM-2013 tables that come with AMCSET. A perfect crystal for every ion.")
    r.say()
    for part in a.only.split(","):
        if part == "dip":
            dip(r)
        elif part == "boron":
            boron(r, [float(x) for x in a.fractions.split(",")])
        else:
            sys.exit(f"unknown part '{part}' (dip, boron)")
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "results.txt").write_text("\n".join(r.lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

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

"""AMCSET's electronic stopping against measurements and reference tables.

  iaea   every measurement of the IAEA electronic stopping power database
         (https://nds.iaea.org/stopping/, the StoppingPower.csv of its full
         download) in an elemental target, against `amcset stopping` with the
         bundled SRIM-2013 tables (srim) and without them (builtin); Ni ions
         also with physics.stopping_corrections.nickel_ions
  nist   NIST PSTAR and ASTAR (protons and alpha particles; ICRU 49) in Si, Fe,
         Cu, W, Au and water: electronic stopping of the srim, builtin and
         bethe-bloch models, and the CSDA range

    python tools/validation/stopping_reference.py --exe build/Release/amcset.exe \\
        --iaea StoppingPower.csv --nist DIR --out build/validation/stopping

DIR holds pstar_<material>.txt and astar_<material>.txt (the NIST text tables:
MeV; electronic, nuclear, total stopping in MeV cm2/g; CSDA and projected range
in g/cm2). The data files are not in the repository (docs/VALIDATION.md, "How the
cases are rerun", says where they come from). Standard library only.
"""

import argparse
import concurrent.futures
import csv
import math
import statistics
import subprocess
import sys
from pathlib import Path

SYMBOLS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
           "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb "
           "Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U").split()
Z = {s: i + 1 for i, s in enumerate(SYMBOLS)}
# IAEA target names of elemental targets that are not plain symbols; diatomic gases.
TARGET_ALIAS = {"Carbon": "C", "H2": "H", "N2": "N", "O2": "O", "F2": "F", "Cl2": "Cl"}
UNITS = "MeV/(mg/cm2)"
# References of the IAEA database that are calculations, not measurements: tables
# (SRIM-2013, MSTAR, ICRU 49), models (DPASS, CasP, ESPNN, SLPA-LM, MELF-GOS, TILDA-V)
# and TDDFT. The SRIM-2013 rows check the bundled tables instead.
CALCULATED = {"SRIM2013", "MSTAR2003", "ESPNN", "DPASS21.06", "ICRU49", "SLPA-LM", "CasP6.0", "TILDA-V(g)"}
CALCULATED_PREFIXES = ("TDDFT", "rtTDDFT", "MELFGOS")


def calculated(ref):
    return ref in CALCULATED or ref.startswith(CALCULATED_PREFIXES)

# --- running AMCSET ------------------------------------------------------------------


def table(exe, ion, target, e_from_keV, e_to_keV, model="srim", gas=False, fix=False, density=None):
    """`amcset stopping` on 50 points per decade: lists (E keV, S_e, S_n, CSDA range um)."""
    cmd = [str(exe), "stopping", "--ion", ion, "--target", target, "--from", f"{e_from_keV:.6g} keV",
           "--to", f"{e_to_keV:.6g} keV", "--per-decade", "50", "--units", UNITS, "--csv", "--model", model]
    if density:
        cmd += ["--density", f"{density} g/cm3"]
    if gas:
        cmd.append("--gas")
    if fix:
        cmd.append("--correct-nickel")
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}: {r.stderr.strip()}")
    rows = []
    for line in r.stdout.splitlines():
        if not line or line.startswith("#"):
            continue
        e, se, sn, _, rng = (float(v) for v in line.split(","))
        rows.append((e, se, sn, rng))
    return rows


def interpolate(rows, e, column):
    """Log-log interpolation of column (1: S_e, 2: S_n, 3: range) at energy e."""
    lo = 0
    hi = len(rows) - 1
    if e <= rows[0][0]:
        return rows[0][column]
    if e >= rows[hi][0]:
        return rows[hi][column]
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if rows[mid][0] <= e:
            lo = mid
        else:
            hi = mid
    (e0, y0), (e1, y1) = (rows[lo][0], rows[lo][column]), (rows[hi][0], rows[hi][column])
    if y0 <= 0 or y1 <= 0:
        return y0 + (y1 - y0) * (e - e0) / (e1 - e0)
    t = math.log(e / e0) / math.log(e1 / e0)
    return math.exp(math.log(y0) + t * (math.log(y1) - math.log(y0)))


def default_mass(exe, symbol, cache={}):
    """Mass of the most abundant isotope as `amcset stopping` uses it (from its header)."""
    if symbol not in cache:
        rows = subprocess.run([str(exe), "stopping", "--ion", symbol, "--target", "Si", "--from", "1 keV",
                               "--to", "10 keV", "--per-decade", "1"], capture_output=True, text=True).stdout
        head = rows.splitlines()[0]
        cache[symbol] = float(head.split("(")[1].split(" amu")[0])
    return cache[symbol]


# --- IAEA database -------------------------------------------------------------------


def read_iaea(path):
    """Elemental targets: dicts (ion, target, gas) -> list of (keV/u, S_e, error %, ref) for the
    measurements and for the database's own copy of SRIM-2013; the number of rows left out
    (compounds, unusable energies) and of calculated rows other than SRIM-2013's."""
    pairs, srim = {}, {}
    skipped = other_calc = 0
    with open(path, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            ion = row["projectile_name"].strip()
            target = TARGET_ALIAS.get(row["target_name"].strip(), row["target_name"].strip())
            phase = row["phase"].strip()
            if ion not in Z or target not in Z or row["stopping_unit_converted"] != UNITS:
                skipped += 1
                continue
            try:
                e = float(row["energy"])
                s = float(row["stopping_power_converted"])
                mass = float(row["ion_isotope"]) if row["ion_isotope"].strip() else None
            except ValueError:
                skipped += 1
                continue
            unit = row["energy_unit"].strip().lower()
            if unit in ("mev/u", "kev/u"):
                e_u = e * (1000.0 if unit == "mev/u" else 1.0)
            elif unit in ("mev", "kev") and mass:
                e_u = e * (1000.0 if unit == "mev" else 1.0) / mass
            else:
                skipped += 1
                continue
            if s <= 0 or e_u <= 0:
                skipped += 1
                continue
            try:
                err = float(row["percentage_error"])
            except ValueError:  # blank, or a range such as "3.9-4.1"
                err = float("nan")
            ref = row["ref_id"].strip()
            if ref == "SRIM2013":
                dest = srim
            elif calculated(ref):
                other_calc += 1
                continue
            else:
                dest = pairs
            dest.setdefault((ion, target, phase == "g"), []).append((e_u, s, err, ref))
    return pairs, srim, skipped, other_calc


BANDS = [("< 0.025", 0.0, 25.0), ("0.025-0.25", 25.0, 250.0), ("0.25-2.5", 250.0, 2500.0),
         ("> 2.5 MeV/u", 2500.0, float("inf"))]


def summary_line(label, ratios):
    if not ratios:
        return f"   {label:28s} {0:6d}"
    med = statistics.median(ratios)
    dev = sorted(abs(r - 1.0) for r in ratios)
    w10 = sum(1 for d in dev if d <= 0.10) / len(dev)
    w20 = sum(1 for d in dev if d <= 0.20) / len(dev)
    return (f"   {label:28s} {len(ratios):6d} {med:8.3f} {dev[len(dev) // 2]:8.3f} {100 * w10:7.1f} %"
            f" {100 * w20:7.1f} %")


def iaea(exe, path, lines, workers):
    pairs, srim_rows, skipped, other_calc = read_iaea(path)
    masses = {ion: default_mass(exe, ion) for ion in {k[0] for k in pairs} | {k[0] for k in srim_rows}}

    jobs = []
    for (ion, target, gas), pts in pairs.items():
        lo = min(p[0] for p in pts) * masses[ion] / 1.2
        hi = max(p[0] for p in pts) * masses[ion] * 1.2
        lo, hi = max(lo, 1e-3), min(hi, 1e8)
        for model, fix in (("srim", False), ("builtin", False)) + ((("srim", True),) if ion == "Ni" else ()):
            jobs.append(((ion, target, gas, model, fix), (exe, ion, target, lo, hi, model, gas, fix)))
    for (ion, target, gas), pts in srim_rows.items():
        if (ion, target, gas) in pairs:
            continue  # the srim table of the measured pair covers it (range widened below)
        lo = max(min(p[0] for p in pts) * masses[ion] / 1.2, 1e-3)
        hi = min(max(p[0] for p in pts) * masses[ion] * 1.2, 1e8)
        jobs.append(((ion, target, gas, "srim", False), (exe, ion, target, lo, hi, "srim", gas, False)))
    # a measured pair's srim table must also span the SRIM-2013 rows of that pair
    for i, (key, args) in enumerate(jobs):
        ion, target, gas, model, fix = key
        if model == "srim" and not fix and (ion, target, gas) in pairs and (ion, target, gas) in srim_rows:
            pts = srim_rows[(ion, target, gas)]
            lo = min(args[3], max(min(p[0] for p in pts) * masses[ion] / 1.2, 1e-3))
            hi = max(args[4], min(max(p[0] for p in pts) * masses[ion] * 1.2, 1e8))
            jobs[i] = (key, (exe, ion, target, lo, hi, model, gas, fix))
    tables = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(table, *args): key for key, args in jobs}
        for fut in concurrent.futures.as_completed(futures):
            key = futures[fut]
            try:
                tables[key] = fut.result()
            except RuntimeError as e:
                print(e, file=sys.stderr)

    # ratio AMCSET / measured for every point
    points = []  # (ion, target, gas, keV/u, measured, srim, builtin, srim_fix, ref)
    for (ion, target, gas), pts in pairs.items():
        ts = tables.get((ion, target, gas, "srim", False))
        tb = tables.get((ion, target, gas, "builtin", False))
        tf = tables.get((ion, target, gas, "srim", True))
        if not ts or not tb:
            continue
        for e_u, s, _, ref in pts:
            e = e_u * masses[ion]
            points.append((ion, target, gas, e_u, s, interpolate(ts, e, 1), interpolate(tb, e, 1),
                           interpolate(tf, e, 1) if tf else None, ref))

    # the database's own copy of SRIM-2013 against the bundled tables
    copy = []
    for (ion, target, gas), pts in srim_rows.items():
        ts = tables.get((ion, target, gas, "srim", False))
        if ts:
            copy += [interpolate(ts, e_u * masses[ion], 1) / s for e_u, s, _, _ in pts]
    dev = sorted(abs(r - 1.0) for r in copy)

    head = f"   {'':28s} {'points':>6s} {'median':>8s} {'med|r-1|':>8s} {'<=10%':>9s} {'<=20%':>9s}"
    lines += ["== IAEA electronic stopping power database (StoppingPower.csv): elemental targets",
              f"   Its SRIM-2013 rows (ref SRIM2013, {len(copy)} values of {len(srim_rows)} pairs) against the "
              f"bundled tables: median ratio {statistics.median(copy):.4f},",
              f"   median |r-1| {dev[len(dev) // 2]:.4f}, 99th percentile {dev[int(0.99 * len(dev))]:.4f}, "
              f"largest {dev[-1]:.4f}.",
              f"   Measurements: {len(points)} of {len(pairs)} ion-target pairs; left out: {other_calc} rows of other "
              "calculations (MSTAR, DPASS, CasP, ESPNN,",
              f"   ICRU 49, TDDFT, ...) and {skipped} rows of compounds, mixtures or without a usable energy.",
              "   AMCSET / measured; srim: the bundled SRIM-2013 tables; builtin: AMCSET's stopping without them.",
              "   Energy bands in MeV/u.",
              head]
    for model, col in (("srim", 5), ("builtin", 6)):
        lines.append(summary_line(f"{model}, all", [p[col] / p[4] for p in points]))
        for label, a, b in BANDS:
            lines.append(summary_line(f"{model}, {label}", [p[col] / p[4] for p in points if a <= p[3] < b]))
    lines.append("")

    # per ion (srim), sorted by Z
    lines += ["   srim tables per ion: median AMCSET / measured in each band (points in brackets)",
              f"   {'ion':4s} " + " ".join(f"{label:>16s}" for label, _, _ in BANDS)]
    for ion in sorted({p[0] for p in points}, key=lambda s: Z[s]):
        cells = []
        for _, a, b in BANDS:
            rs = [p[5] / p[4] for p in points if p[0] == ion and a <= p[3] < b]
            cells.append(f"{statistics.median(rs):7.3f} ({len(rs):5d})" if rs else f"{'':16s}")
        lines.append(f"   {ion:4s} " + " ".join(f"{c:>16s}" for c in cells))
    lines.append("")

    # Ni and its neighbours, by target, in the band of SRIM's Ni dip
    lines += ["   Ni ions and their neighbours, 0.02-0.3 MeV/u (where SRIM-2013's Ni stopping dips): median",
              "   AMCSET / measured by target (points); Ni also with physics.stopping_corrections.nickel_ions",
              f"   {'ion':4s} {'target':7s} {'points':>6s} {'srim':>7s} {'nickel fix':>10s}"]
    for ion in ("Fe", "Co", "Ni", "Cu"):
        for target in sorted({p[1] for p in points if p[0] == ion}, key=lambda s: Z[s]):
            sel = [p for p in points if p[0] == ion and p[1] == target and 20.0 <= p[3] <= 300.0]
            if not sel:
                continue
            m = statistics.median(p[5] / p[4] for p in sel)
            fx = f"{statistics.median(p[7] / p[4] for p in sel):10.3f}" if ion == "Ni" else ""
            lines.append(f"   {ion:4s} {target:7s} {len(sel):6d} {m:7.3f} {fx}")
        sel = [p for p in points if p[0] == ion and 20.0 <= p[3] <= 300.0]
        if sel:
            m = statistics.median(p[5] / p[4] for p in sel)
            fx = f"{statistics.median(p[7] / p[4] for p in sel):10.3f}" if ion == "Ni" else ""
            lines.append(f"   {ion:4s} {'(all)':7s} {len(sel):6d} {m:7.3f} {fx}")
    lines.append("   Ni ions at every energy, by band (srim / with the fix):")
    for label, a, b in BANDS:
        sel = [p for p in points if p[0] == "Ni" and a <= p[3] < b]
        if sel:
            lines.append(f"      {label:12s} {len(sel):5d} {statistics.median(p[5] / p[4] for p in sel):7.3f} "
                         f"{statistics.median(p[7] / p[4] for p in sel):7.3f}")
    lines.append("")
    return points


# --- NIST PSTAR / ASTAR --------------------------------------------------------------

NIST_MATERIALS = [("silicon", "Si", None), ("iron", "Fe", None), ("copper", "Cu", None), ("tungsten", "W", None),
                  ("gold", "Au", None), ("water", "H:2,O:1", "1.0")]


def read_nist(path):
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split()
        if len(parts) < 6:
            continue
        try:
            vals = [float(v) for v in parts[:6]]
        except ValueError:
            continue
        rows.append(vals)  # E MeV, S_e, S_n, S_tot (MeV cm2/g), CSDA (g/cm2), projected (g/cm2)
    return rows


def nist(exe, folder, lines):
    lines += ["== NIST PSTAR / ASTAR (ICRU 49): AMCSET electronic stopping / NIST, and CSDA range / NIST",
              "   at the NIST energies; srim: bundled SRIM-2013 tables, builtin: without them, bb: bethe-bloch",
              "   (only above 1 MeV/u). Columns: median ratio and largest deviation over the energies named."]
    for prog, ion, mass in (("pstar", "H", 1.007), ("astar", "He", 4.0026)):
        lines.append(f"   {prog.upper()} ({ion} ions)")
        lines.append(f"   {'material':9s} {'band (MeV)':>14s} {'n':>4s} {'srim':>16s} {'builtin':>16s} {'bb':>16s}"
                     f" {'CSDA srim':>16s}")
        for name, target, density in NIST_MATERIALS:
            f = folder / f"{prog}_{name}.txt"
            if not f.exists():
                lines.append(f"   {name:9s} (no file)")
                continue
            ref = read_nist(f)
            lo, hi = ref[0][0] * 1000.0 / 1.2, ref[-1][0] * 1000.0 * 1.2
            t = {m: table(exe, ion, target, lo, hi, m, density=density) for m in ("srim", "builtin", "bethe-bloch")}
            # mass stopping: amcset's MeV/(mg/cm2) = 1000 MeV cm2/g
            for band, a, b in (("0.01-0.1", 0.01, 0.1), ("0.1-1", 0.1, 1.0), ("1-10", 1.0, 10.0),
                               ("10-100", 10.0, 100.0), ("100-1000", 100.0, 1000.0)):
                sel = [r for r in ref if a <= r[0] / mass < b] if ion == "He" else [r for r in ref if a <= r[0] < b]
                if not sel:
                    continue
                cells = []
                for m in ("srim", "builtin", "bethe-bloch"):
                    if m == "bethe-bloch" and a < 1.0:
                        cells.append(f"{'':16s}")
                        continue
                    rs = [interpolate(t[m], r[0] * 1000.0, 1) * 1000.0 / r[1] for r in sel]
                    worst = max(rs, key=lambda x: abs(x - 1.0))
                    cells.append(f"{statistics.median(rs):7.3f} ({worst:6.3f})")
                rho = float(density) if density else densities[target]
                rs = [interpolate(t["srim"], r[0] * 1000.0, 3) / (r[4] / rho * 1e4) for r in sel]
                worst = max(rs, key=lambda x: abs(x - 1.0))
                cells.append(f"{statistics.median(rs):7.3f} ({worst:6.3f})")
                lines.append(f"   {name:9s} {band:>14s} {len(sel):4d} " + " ".join(cells))
        lines.append("")
    lines.append("   (He bands are in MeV/u; in brackets the ratio farthest from 1 in the band; water: H2O at")
    lines.append("   1 g/cm3 by Bragg's rule, which SRIM's tables do not correct for the phase or the bonds.)")
    lines.append("")


densities = {}


def read_densities(exe):
    """Densities amcset stopping uses for the NIST elements (its header), g/cm3."""
    for _, target, d in NIST_MATERIALS:
        if d is not None:
            continue
        head = subprocess.run([str(exe), "stopping", "--ion", "H", "--target", target, "--from", "1 MeV",
                               "--to", "2 MeV", "--per-decade", "1"], capture_output=True, text=True).stdout
        densities[target] = float(head.split("density ")[1].split(" g/cm3")[0])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--iaea", type=Path, help="StoppingPower.csv of the IAEA database")
    ap.add_argument("--nist", type=Path, help="folder with pstar_*.txt and astar_*.txt")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()
    lines = ["AMCSET's electronic stopping against measurements and reference tables "
             "(tools/validation/stopping_reference.py)", ""]
    a.out.mkdir(parents=True, exist_ok=True)
    if a.iaea:
        points = iaea(a.exe, a.iaea, lines, a.workers)
        with open(a.out / "iaea_ratios.csv", "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f)
            w.writerow(["ion", "target", "gas", "keV_per_u", "measured", "srim", "builtin", "srim_nickel_fix", "ref"])
            for p in points:
                w.writerow(p)
    if a.nist:
        read_densities(a.exe)
        nist(a.exe, a.nist, lines)
    (a.out / "results.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

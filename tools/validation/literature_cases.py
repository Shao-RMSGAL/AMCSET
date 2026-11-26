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

"""AMCSET (SRIM mode, amorphous targets) against measurements from the literature.

Each case states its source, the measured values and the conditions next to
the scenario that reproduces them; docs/VALIDATION.md quotes results.txt.

  ni_ranges     Ni ions in Fe and Mo: SRIM-2013 overestimates their ranges
                (Gurbich et al. 2025); AMCSET with and without
                physics.stopping_corrections.nickel_ions
  ni_peaks      measured peak depths of implanted Ni in Ni, Fe and two steels, 2-11.5 MeV
                (Chauhan et al. 2026, Pena et al. 2025, Song et al. 2026), with and without
                the correction, and the shift by sputtering during the implantation
  ranges        mean projected ranges and straggling (or peaks) of range_data.py: Behar et
                al. 1995 (Si, SiO2, C, B, SiC; H, Li, B in Si), Oosterhoff 1988 (B, P in Si),
                Zhang et al. 2009 (Pb, Au in SiC); OUT/ranges.csv
  sputtering    sputtering yields of sputter_data.py (Ar on Si, Fe, Cu, Ag, W, Au; H, D, He
                on Fe and W; Yamamura et al.); OUT/sputtering.csv
  published_srim  published SRIM damage calculations (Stoller et al. 2013, Weber and Zhang
                2019): vacancies and damage energy per ion
  as_si_channeling  the crystal transport: channeling tails of 15 keV As in (100) Si against
                tilt (Tasch et al. 1994)
  sb_si         50 keV Sb into Si with its native oxide, SIMS (Titze et al. 2021)
  quenching     ionization yield of Ge recoils in Ge (Bonhomme et al. 2022) and of
                Si recoils in Si (Chavarria et al. 2016, Izraelevitch et al. 2017),
                with Lindhard's prediction

    python tools/validation/literature_cases.py --exe build/Release/amcset.exe \
        --out build/validation/literature [--only ni_ranges,ni_peaks,...] [--seeds 4]

Runs are cached in OUT/runs (delete a folder to rerun it). Standard library only;
about half an hour for every case on the workstation.
"""

import argparse
import json
import math
import random
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))  # range_data, sputter_data (also under python -I)

# --- measured data (with sources) -------------------------------------------------

# A.F. Gurbich, P.S. Prusachenko, T.L. Bobrovskiy, M.V. Bokhovko, "On the SRIM
# prediction of ranges for Ni ions", Appl. Radiat. Isot. 215, 111590 (2025),
# doi:10.1016/j.apradiso.2024.111590, abstract: SRIM-2013 overestimates the
# projected range of 5 MeV Ni in Mo by 1.3 times and of 11.5 MeV Ni in Fe by
# 1.5 times; TRIM-98 with radiation-enhanced diffusion agrees with the data.
NI_RANGES = [
    # (target, density g/cm3, Z2, energy MeV, SRIM-2013 / measured)
    ("Mo", 10.22, 42, 5.0, 1.3),
    ("Fe", 7.8658, 26, 11.5, 1.5),
]

# Measured peak depths (um) of implanted Ni, against the SRIM-2013 peaks the papers give.
#  C: R. Chauhan et al., "A study of electronic-stopping deficiencies in SRIM using
#     high-energy nickel ion implantation into pure nickel and pure iron", J. Nucl. Mater. 633,
#     157044 (2026), doi:10.1016/j.jnucmat.2026.157044: 62Ni, normal incidence, room
#     temperature, ToF-SIMS (bi-Gaussian peak); SRIM-2013 quick, Ed 40 eV. Fe values from the
#     table of S below (which quotes C); 5 MeV in Ni from C's text: the SRIM peak is 625 nm
#     (about 36 %) deeper than the measured one.
#  P: M. Pena et al., Nucl. Instrum. Methods B 565, 165738 (2025),
#     doi:10.1016/j.nimb.2025.165738: 62Ni into Ni, 4.6e17 /cm2, ToF-SIMS peaks.
#  S: M. Song, G.S. Was, F.A. Garner, J. Nucl. Mater. 633, 157002 (2026),
#     doi:10.1016/j.jnucmat.2026.157002, Table 2 (peak depths; their SRIM-2013 bins 40-50 nm):
#     pure Fe (C's SIMS; 8 MeV: Lin et al., EDX in Fe and Fe-10Cr; 11.5 MeV: Gurbich et al.
#     2025, resonance profiling), 8 MeV in 316 SS and 9 MeV in 304 SS (STEM-EDX).
# The measured peaks are not corrected for sputtering during the implantation (C, §2.2).
# (key, ion mass amu, target, energy MeV, measured um, uncertainty um, paper's SRIM um,
#  fluence /cm2 or None, source)
NI_PEAKS = [
    ("Ni_2.0", 61.928, "Ni", 2.0, 0.364, 0.05, 0.720, 4.6e17, "P"),
    ("Ni_3.5", 61.928, "Ni", 3.5, 1.008, 0.05, 1.290, 4.6e17, "P"),
    ("Ni_5.0", 61.928, "Ni", 5.0, 1.111, 0.05, 1.736, 1.5e17, "C"),
    ("Fe_2.0", 61.928, "Fe", 2.0, 0.584, 0.010, 0.785, 2.32e17, "C/S"),
    ("Fe_3.5", 61.928, "Fe", 3.5, 1.137, 0.010, 1.380, 4.65e17, "C/S"),
    ("Fe_5.0", 61.928, "Fe", 5.0, 1.389, 0.010, 1.890, 4.64e17, "C/S"),
    ("Fe_8.0", 57.935, "Fe", 8.0, 1.750, None, 2.50, None, "S (Lin)"),
    ("Fe_11.5", 57.935, "Fe", 11.5, 1.92, 0.16, 2.80, None, "S (Gurbich)"),
    ("316_8.0", 57.935, "316", 8.0, 1.815, 0.003, 2.44, 1.22e17, "S"),
    ("304_9.0", 57.935, "304", 9.0, 1.925, 0.030, 2.62, 1.29e17, "S"),
]
# Targets: density (g/cm3) and composition (at.%; S Table 1, the minor elements C, P, S, N,
# Co, Cu left out). The steels' densities are not in S: handbook values.
NI_TARGETS = {
    "Ni": (8.902, [("Ni", 100.0)]),
    "Fe": (7.8658, [("Fe", 100.0)]),
    "316": (7.99, [("Fe", 64.53), ("Cr", 18.7), ("Ni", 12.63), ("Mn", 1.73), ("Mo", 1.56), ("Si", 0.85)]),
    "304": (7.93, [("Fe", 65.71), ("Cr", 21.01), ("Ni", 10.08), ("Mn", 1.81), ("Mo", 0.30), ("Si", 1.09)]),
}

# M. Titze et al., "Evaluation of the accuracy of stopping and range of ions in
# matter simulations through secondary ion mass spectrometry and Rutherford
# backscattering spectrometry for low energy heavy ion implantation", J. Vac.
# Sci. Technol. A 39, 063222 (2021), doi:10.1116/6.0001406 (OSTI 1830531,
# SAND2021-13560J). 50 keV Sb, 7 deg, Si with a 2.9 +- 1.2 nm native oxide
# (SiO2 2.32 g/cm3, Si 2.329 g/cm3; Ed/El/Es O 28/3/2, Si 15/2/4.7 eV,
# full cascades, 99999 ions): SIMS depth 27.5 +- 1.3 nm, straggle 9.7 nm;
# their SRIM run 32.8 nm.
SB_SI = {"energy_keV": 50.0, "angle_deg": 7.0, "oxide_nm": 2.9, "measured_nm": 27.5, "measured_err_nm": 1.3,
         "straggle_nm": 9.7, "srim_nm": 32.8}

# Ionization yield (ionization energy / recoil energy) of nuclear recoils.
# Ge: A. Bonhomme et al., Eur. Phys. J. C 82, 815 (2022), arXiv:2202.03754,
# Table 5 (E_nr keV, q, total uncertainty); compatible with Lindhard k = 0.162.
GE_YIELD = [(0.97, 0.181, 0.013), (1.17, 0.192, 0.019), (1.54, 0.190, 0.011), (1.58, 0.193, 0.010),
            (1.92, 0.190, 0.020), (1.99, 0.187, 0.006)]
# Si: A.E. Chavarria et al., Phys. Rev. D 94, 082007 (2016), arXiv:1608.00957,
# Table 1 (E_e keVee, E_r keVnr, systematic uncertainty of E_r).
SI_CHAVARRIA = [(0.06, 0.68, 0.10), (0.12, 1.05, 0.10), (0.18, 1.36, 0.08), (0.24, 1.66, 0.06),
                (0.30, 1.94, 0.06), (0.39, 2.28, 0.03)]
# Si: F. Izraelevitch et al., JINST 12, P06014 (2017), arXiv:1702.00873,
# Table 1 (E_NR keVnr, efficiency, systematic uncertainty).
SI_IZRAELEVITCH = [(1.79, 0.168, 0.034), (2.97, 0.215, 0.020), (4.87, 0.258, 0.011), (7.67, 0.280, 0.006),
                   (10.77, 0.308, 0.005), (14.48, 0.329, 0.004), (20.67, 0.360, 0.004)]


def lindhard_yield(e_keV, z, k):
    """Lindhard's ionization yield k g / (1 + k g), g = 3 eps^0.15 + 0.7 eps^0.6 + eps,
    eps = 11.5 E/keV Z^(-7/3) (Lindhard et al., Mat. Fys. Medd. Dan. Vid. Selsk. 33, 10
    (1963), in the form of Lewin and Smith, Astropart. Phys. 6, 87 (1996))."""
    eps = 11.5 * e_keV * z ** (-7.0 / 3.0)
    g = 3.0 * eps ** 0.15 + 0.7 * eps ** 0.6 + eps
    return k * g / (1.0 + k * g)


# --- running AMCSET ----------------------------------------------------------------


def run(exe, scenario, work, extra=()):
    """Run a scenario (cached in work/), return its summary.json."""
    work.mkdir(parents=True, exist_ok=True)
    summary = work / "out" / "summary.json"
    if not summary.exists():
        scenario = dict(scenario)
        scenario["output"] = dict(scenario.get("output", {}), directory=str((work / "out").resolve()),
                                  files={"formats": ["trim"]})
        sc = work / "scenario.json"
        sc.write_text(json.dumps(scenario, indent=1), encoding="utf-8")
        r = subprocess.run([str(exe), "run", str(sc), "-q", *extra], capture_output=True, text=True)
        if r.returncode != 0 or not summary.exists():
            sys.exit(f"{work}: amcset failed\n{r.stdout}\n{r.stderr}")
    return json.loads(summary.read_text(encoding="utf-8"))


# Structured results of every case (OUT/results.json, read by make_figures.py).
RESULTS = {}


def record(case, **values):
    RESULTS.setdefault(case, []).append(values)


def mean_sem(values):
    n = len(values)
    m = sum(values) / n
    if n < 2:
        return m, float("nan")
    return m, math.sqrt(sum((v - m) ** 2 for v in values) / (n - 1) / n)


def element(sym, ed=None, eb=None, es=None, fraction=1.0):
    e = {"element": sym, "stoichiometry": fraction}
    if ed is not None:
        e["Ed"] = f"{ed} eV"
    if eb is not None:
        e["Eb"] = f"{eb} eV"
    if es is not None:
        e["Es"] = f"{es} eV"
    return e


def ni_ranges(exe, runs, seeds, lines):
    lines += ["== Ni ions in Mo and Fe (Gurbich et al. 2025): mean projected range",
              "   quick (Kinchin-Pease) calculation, 2000 ions x seeds, bundled SRIM-2013 tables; the measured",
              "   range is SRIM-2013's divided by the paper's factor.",
              f"   {'case':22s} {'SRIM-2013 tables':>18s} {'nickel_ions fix':>18s} {'ratio':>7s}  {'paper: SRIM/measured':>20s}"]
    for target, rho, _, e_MeV, factor in NI_RANGES:
        res = {}
        for fix in (False, True):
            vals = []
            for k in range(seeds):
                sc = {"mode": "srim", "title": f"Ni {e_MeV} MeV in {target}",
                      "ion": {"element": "Ni", "energy": f"{e_MeV} MeV"},
                      "layers": [{"name": target, "thickness": "8 um", "density": f"{rho} g/cm3",
                                  "elements": [element(target)]}],
                      "ions": 2000, "seed": 101 + k, "threads": 0, "damage": "quick",
                      "physics": {"stopping_corrections": {"nickel_ions": fix}}}
                s = run(exe, sc, runs / f"ni_{target}_{e_MeV:g}MeV_{'fix' if fix else 'srim'}_s{k}")
                vals.append(s["range"]["mean_A"] / 1e4)
            res[fix] = mean_sem(vals)
        (a, ea), (b, eb_) = res[False], res[True]
        record("ni_ranges", target=target, energy_MeV=e_MeV, srim_um=a, srim_sem=ea, fix_um=b, fix_sem=eb_,
               paper_factor=factor)
        lines.append(f"   {e_MeV:g} MeV Ni -> {target:4s}       {a:8.3f} +- {ea:.3f} um {b:8.3f} +- {eb_:.3f} um"
                     f" {a / b:7.3f}  {factor:20.2f}")
    lines.append("")


def ion_profile(work):
    """Stopped-ion histogram of a run (RANGE.txt): bin centres (A), counts per bin (any scale),
    and the target's atom density (atoms/cm3, of the first layer)."""
    s = json.loads((work / "out" / "summary.json").read_text(encoding="utf-8"))
    b = s["depth_bins"]
    text = (work / "out" / "RANGE.txt").read_text(encoding="utf-8", errors="replace").splitlines()
    density = None
    counts = []
    started = False
    for line in text:
        if density is None and "Density =" in line and "atoms/cm3" in line:
            density = float(line.split("Density =")[1].split("atoms/cm3")[0])
        if line.startswith("-----------  ----------"):  # the rule under the table's column heads
            started = True
            continue
        if started:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    counts.append(float(parts[1]))
                except ValueError:
                    break
    centres = [b["start_A"] + (i + 0.5) * b["width_A"] for i in range(len(counts))]
    return centres, counts, density


def peak_depth(centres, counts):
    """Depth of the profile's maximum: a parabola through ln(y) of the 1-2-3-2-1 smoothed profile
    where it exceeds 60 % of its maximum (a Gaussian's peak; the papers fit Gaussians or splines)."""
    w = [1, 2, 3, 2, 1]
    n = len(counts)
    sm = []
    for i in range(n):
        acc = norm = 0.0
        for k, wk in enumerate(w):
            j = i + k - 2
            if 0 <= j < n:
                acc += wk * counts[j]
                norm += wk
        sm.append(acc / norm)
    i0 = max(range(n), key=lambda i: sm[i])
    lo = hi = i0
    while lo > 0 and sm[lo - 1] > 0.6 * sm[i0]:
        lo -= 1
    while hi < n - 1 and sm[hi + 1] > 0.6 * sm[i0]:
        hi += 1
    if hi - lo < 2:
        return centres[i0]
    # least-squares parabola y = a + b x + c x^2 on ln(sm)
    xs = [centres[i] for i in range(lo, hi + 1)]
    ys = [math.log(sm[i]) for i in range(lo, hi + 1)]
    x0 = centres[i0]
    sx = [sum((x - x0) ** k for x in xs) for k in range(5)]
    sy = [sum(y * (x - x0) ** k for x, y in zip(xs, ys)) for k in range(3)]
    m = [[sx[0], sx[1], sx[2]], [sx[1], sx[2], sx[3]], [sx[2], sx[3], sx[4]]]

    def det(a):
        return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
                + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))

    d = det(m)
    coef = []
    for col in range(3):
        mc = [row[:] for row in m]
        for r in range(3):
            mc[r][col] = sy[r]
        coef.append(det(mc) / d)
    _, b1, c2 = coef
    if c2 >= 0:
        return centres[i0]
    return x0 - b1 / (2.0 * c2)


def eroded(centres, counts, d):
    """The implanted profile seen from the final surface when a depth d (A) was sputtered away at a
    constant rate during the implantation: the mean of the range profile shifted by 0 to d."""
    if d <= 0:
        return counts
    width = centres[1] - centres[0]

    def p(x):
        t = (x - centres[0]) / width
        i = int(math.floor(t))
        if i < 0 or i + 1 >= len(counts):
            return 0.0
        return counts[i] + (counts[i + 1] - counts[i]) * (t - i)

    steps = 64
    return [sum(p(x + d * (k + 0.5) / steps) for k in range(steps)) / steps for x in centres]


def ni_peaks(exe, runs, seeds, lines):
    lines += ["== Peak depth of implanted Ni (Chauhan et al. 2026, Pena et al. 2025, Song et al. 2026)",
              "   quick (Kinchin-Pease), Ed 40 eV, 20000 ions x seeds, bundled SRIM-2013 tables; 'fix':",
              "   physics.stopping_corrections.nickel_ions. Sputtering: full cascades, 1000 ions x 2 (with the fix);",
              "   d = Y x fluence / N is the depth sputtered away during the implantation, and 'fix, sputtered'",
              "   the peak of the fixed profile seen from the final surface (its shifts by 0 to d averaged).",
              f"   {'case':9s} {'src':11s} {'measured':>13s} {'paper SRIM':>10s} {'AMCSET':>15s} {'/paper':>6s}"
              f" {'fix':>15s} {'/meas':>6s} {'Y':>5s} {'d nm':>5s} {'fix, sputt.':>11s} {'/meas':>6s}"]
    for key, mass, target, e_MeV, meas, err, srim_paper, fluence, src in NI_PEAKS:
        rho, comp = NI_TARGETS[target]
        thick_um = 2.2 * srim_paper

        def scenario(fix, k, ions, damage):
            return {"mode": "srim", "title": f"Ni {e_MeV} MeV in {target}",
                    "ion": {"element": "Ni", "mass_amu": mass, "energy": f"{e_MeV} MeV"},
                    "layers": [{"name": target, "thickness": f"{thick_um:.4g} um", "density": f"{rho} g/cm3",
                                "elements": [element(sym, 60 if sym == "Mo" else 40, None, None, f)
                                             for sym, f in comp]}],
                    "ions": ions, "seed": 101 + k, "threads": 0, "damage": damage,
                    "physics": {"stopping_corrections": {"nickel_ions": fix}}}

        res = {}
        hist = {}
        for fix in (False, True):
            peaks = []
            total = None
            for k in range(seeds):
                work = runs / f"nipeak_{key}_{'fix' if fix else 'srim'}_s{k}"
                run(exe, scenario(fix, k, 20000, "quick"), work)
                c, n, density = ion_profile(work)
                peaks.append(peak_depth(c, n) / 1e4)
                total = n if total is None else [a + b for a, b in zip(total, n)]
            res[fix] = mean_sem(peaks)
            hist[fix] = (c, total, density)
        y = []
        for k in range(2):
            s = run(exe, scenario(True, k, 1000, "full"), runs / f"nipeak_{key}_sputter_s{k}")
            y.append(sum(e["sputtering_yield"] for e in s["elements"]))
        ym, _ = mean_sem(y)
        c, total, density = hist[True]
        pe = d_A = None
        if fluence:
            d_A = ym * fluence / density * 1e8
            pe = peak_depth(c, eroded(c, total, d_A)) / 1e4
            sput = f"{ym:5.2f} {d_A / 10.0:5.0f} {pe:11.3f} {pe / meas:6.3f}"
        else:
            sput = f"{ym:5.2f} {'-':>5s} {'-':>11s} {'-':>6s}"
        (a, ea), (b, eb_) = res[False], res[True]
        record("ni_peaks", key=key, mass_amu=mass, target=target, energy_MeV=e_MeV, measured_um=meas,
               measured_err_um=err, paper_srim_um=srim_paper, fluence_cm2=fluence, source=src,
               srim_um=a, srim_sem=ea, fix_um=b, fix_sem=eb_, sputtering_yield=ym,
               sputtered_nm=d_A / 10.0 if d_A else None, fix_sputtered_um=pe,
               profile_depth_um=[x / 1e4 for x in c], profile_srim=hist[False][1], profile_fix=total)
        m_txt = f"{meas:.3f} +- {err:.3f}" if err else f"{meas:.3f}       "
        lines.append(f"   {key:9s} {src:11s} {m_txt:>13s} {srim_paper:10.3f} {a:7.3f} +- {ea:.3f} {a / srim_paper:6.3f}"
                     f" {b:7.3f} +- {eb_:.3f} {b / meas:6.3f} {sput}")
    lines += ["   (measured: peak depth in um with the paper's uncertainty; Ni_5.0 from C's text, 625 nm",
              "   shallower than their SRIM peak of about 1.74 um; 8 MeV in Fe and Fe-10Cr has no uncertainty)", ""]


def ranges(exe, runs, seeds, lines):
    """Measured ranges of range_data.GROUPS; every point also in OUT/ranges.csv."""
    import range_data
    seeds = min(seeds, 2)
    rows = []
    lines += ["== Ranges against measurements (tools/validation/range_data.py: sources, targets, methods)",
              f"   quick (Kinchin-Pease), bundled SRIM-2013 tables, 4000 ions x {seeds} (peaks: 20000 x {seeds});",
              "   ratio AMCSET / measured of the mean projected range R_p and of its straggling, or of the peak",
              f"   {'group':11s} {'points':>6s} {'R_p: median':>11s} {'mean':>6s} {'sd':>6s} {'<=10%':>6s}"
              f" {'straggling':>10s}  per ion: median R_p ratio"]
    for key, title, (rho, comp), angle, quantity, points in range_data.GROUPS:
        res = []
        for pt in points:
            ion, mass, e_keV, meas, strag = pt[:5]
            ang = pt[5] if angle is None else angle
            thick_nm = 6.0 * meas + 20.0
            vals, svals = [], []
            for k in range(seeds):
                sc = {"mode": "srim", "title": f"{ion} {e_keV} keV in {key}",
                      "ion": {"element": ion, "mass_amu": mass, "energy": f"{e_keV} keV", "angle": f"{ang} deg"},
                      "layers": [{"name": key, "thickness": f"{thick_nm:.5g} nm", "density": f"{rho} g/cm3",
                                  "elements": [element(sym, fraction=f) for sym, f in comp]}],
                      "ions": 20000 if quantity == "peak" else 4000, "seed": 101 + k, "threads": 0,
                      "damage": "quick"}
                work = runs / f"range_{key}_{ion}{mass:.0f}_{e_keV:g}keV_s{k}"
                s = run(exe, sc, work)
                if quantity == "peak":
                    c, n, _ = ion_profile(work)
                    vals.append(peak_depth(c, n) / 10.0)
                else:
                    vals.append(s["range"]["mean_A"] / 10.0)
                    svals.append(s["range"]["straggling_A"] / 10.0)
            v, _ = mean_sem(vals)
            sv = mean_sem(svals)[0] if svals else None
            r = v / meas
            rs = sv / strag if (sv is not None and strag) else None
            res.append((ion, e_keV, r, rs))
            rows.append((key, ion, mass, e_keV, ang, quantity, meas, strag, round(v, 3),
                         round(sv, 3) if sv is not None else "", round(r, 4), round(rs, 4) if rs else ""))
        rr = [x[2] for x in res]
        ss = [x[3] for x in res if x[3] is not None]
        m = sum(rr) / len(rr)
        sd = math.sqrt(sum((x - m) ** 2 for x in rr) / max(1, len(rr) - 1))
        w10 = sum(1 for x in rr if abs(x - 1.0) <= 0.10) / len(rr)
        per_ion = {}
        for ion, _, r, _ in res:
            per_ion.setdefault(ion, []).append(r)
        ions_txt = " ".join(f"{ion} {sorted(v)[len(v) // 2]:.2f}" for ion, v in per_ion.items())
        strag_txt = f"{sorted(ss)[len(ss) // 2]:10.3f}" if ss else f"{'-':>10s}"
        lines.append(f"   {key:11s} {len(rr):6d} {sorted(rr)[len(rr) // 2]:11.3f} {m:6.3f} {sd:6.3f} {100 * w10:5.0f}%"
                     f" {strag_txt}  {ions_txt}")
        lines.append(f"      {title}")
        if key == "zhang_sic":
            import range_data as rd
            comp_txt = " ".join(f"{ion}{e:g}:{v * pt_meas / rd.ZHANG_SRIM_PEAK[(ion, e)]:.2f}"
                                for (ion, e, v, _), pt_meas in zip(res, [p[3] for p in points]))
            lines.append(f"      AMCSET peak / their SRIM-2008 peak: {comp_txt}")
    out = runs.parent / "ranges.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("group,ion,mass_amu,energy_keV,angle_deg,quantity,measured_nm,measured_straggling_nm,"
                "amcset_nm,amcset_straggling_nm,ratio,straggling_ratio\n")
        for row in rows:
            f.write(",".join(str(x) for x in row) + "\n")
    lines += ["   (sd: spread of the ratios; straggling: median ratio of the straggling where measured)", ""]


def sputtering(exe, runs, seeds, lines):
    """Sputtering yields against sputter_data.MEASURED; every point also in OUT/sputtering.csv."""
    import sputter_data
    rows = []
    for ion, target, e_keV, y_meas, fig, flag in sputter_data.MEASURED:
        # enough ions for about 2000 sputtered atoms (at most 200000); one seed: the counting error
        ions = int(min(200000, max(2000, 2000.0 / max(y_meas, 1e-6))))
        ion_spec = {"element": "H", "mass_amu": 2.0141} if ion == "D" else {"element": ion}
        sc = {"mode": "srim", "title": f"{ion} {e_keV} keV on {target}",
              "ion": dict(ion_spec, energy=f"{e_keV} keV"),
              "layers": [{"name": target, "thickness": "1 um", "elements": [element(target)]}],
              "ions": ions, "seed": 101, "threads": 0, "damage": "full",
              "output": {"range_3d": False, "exit_lists": False}}
        s = run(exe, sc, runs / f"sputter_{ion}_{target}_{e_keV:g}keV", ("--auto-depth",))
        y = sum(el["sputtering_yield"] for el in s["elements"])
        n_sp = y * ions
        err = y / math.sqrt(n_sp) if n_sp > 0 else float("nan")
        rows.append((ion, target, e_keV, y_meas, y, err, fig, flag))

    lines += ["== Sputtering yields at normal incidence: AMCSET (full cascades, TRIM's default Ed, Eb, Es, about",
              "   2000 sputtered atoms per point) against measurements read off Yamamura and Tawara's figures",
              "   (NIFS-DATA-23, -31; reading error 20-30 %, more near the threshold)",
              f"   {'ion-target':11s} {'points':>6s} {'AMCSET/measured':>16s} {'E >= 0.5 keV':>13s} {'range':>13s}"
              f" {'no yield at':>14s}"]
    pairs = []
    for r in rows:
        if (r[0], r[1]) not in pairs:
            pairs.append((r[0], r[1]))

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else None

    all_r, all_hi = [], []
    for ion, target in pairs:
        sel = [r for r in rows if r[0] == ion and r[1] == target and not r[7]]
        a = [r[4] / r[3] for r in sel]
        hi = [r[4] / r[3] for r in sel if r[2] >= 0.5]
        zero = [r[2] for r in sel if r[4] == 0.0]
        all_r += a
        all_hi += hi
        rng = f"{min(hi):.2f}-{max(hi):.2f}" if hi else "-"
        z = f"<= {max(zero):g} keV" if zero else "-"
        lines.append(f"   {ion + ' ' + target:11s} {len(sel):6d} {med(a):16.2f} {med(hi):13.2f} {rng:>13s} {z:>14s}")
    within = sum(1 for x in all_hi if 0.5 <= x <= 2.0) / len(all_hi)
    lines.append(f"   {'all':11s} {len(all_r):6d} {med(all_r):16.2f} {med(all_hi):13.2f}  (E >= 0.5 keV: "
                 f"{100 * within:.0f} % of {len(all_hi)} points within a factor 2)")
    lines += ["   (medians over the energies of each pair; range: smallest and largest ratio at E >= 0.5 keV; 'no",
              "   yield at': energies where AMCSET sputters nothing but a yield was measured; the two points flagged",
              "   as low outliers in sputter_data.py are left out; every point is in sputtering.csv)", ""]
    out = runs.parent / "sputtering.csv"
    with open(out, "w", encoding="utf-8", newline="") as f:
        f.write("ion,target,energy_keV,measured,amcset,amcset_counting_error,figure,flag\n")
        for r in rows:
            f.write(",".join(f"{x:.5g}" if isinstance(x, float) else str(x) for x in r) + "\n")


# Published SRIM calculations of damage (not measurements: a check that AMCSET's SRIM mode
# gives what other SRIM users got).
#  Stoller et al., Nucl. Instrum. Methods B 310, 75 (2013), Table 4 with the erratum (NIM B 459,
#  196 (2019)): SRIM-2008, pure Fe, Ed 40 eV, lattice binding 0, vacancies per ion from VACANCY.txt
#  ("summary") and the damage energy T_dam of their Eq. (5b), keV.
#  Weber and Zhang, Curr. Opin. Solid State Mater. Sci. 23, 100757 (2019), Table 1: full-cascade
#  vacancies nu_FC; Ed 40 eV (Fe, Cr, Ni), Si 35, C 20, O 45, Ti 70, Sr 80 eV; lattice binding 0;
#  densities Fe-10Cr 7.799, Ni 8.912, SiC 3.21, SrTiO3 5.12 g/cm3; SRIM version not stated.
# (source, ion, MeV, target key, damage, published vacancies per ion, published T_dam keV or None)
PUBLISHED_SRIM = [
    ("Stoller T4", "Fe", 0.069, "Fe", "quick", 469, 47.171), ("Stoller T4", "Fe", 0.069, "Fe", "full", 953, 50.097),
    ("Stoller T4", "Fe", 5.0, "Fe", "quick", 7893, 793.255), ("Stoller T4", "Fe", 5.0, "Fe", "full", 16333, 875.652),
    ("Stoller T4", "H", 2.0, "Fe", "quick", 14, 2.512), ("Stoller T4", "H", 2.0, "Fe", "full", 27, 2.602),
    ("Weber T1", "Fe", 10.0, "FeCr", "full", 19197, None), ("Weber T1", "Ni", 0.05, "Ni", "full", 877, None),
    ("Weber T1", "Ni", 5.0, "Ni", "full", 28148, None), ("Weber T1", "Au", 3.0, "Ni", "full", 33214, None),
    ("Weber T1", "Au", 3.0, "SiC", "full", 14577, None), ("Weber T1", "Au", 0.9, "SrTiO3", "full", 4339, None),
]
PUBLISHED_TARGETS = {  # density g/cm3, [(element, stoichiometry, Ed eV)]
    "Fe": (7.866, [("Fe", 1, 40)]), "FeCr": (7.799, [("Fe", 0.9, 40), ("Cr", 0.1, 40)]),
    "Ni": (8.912, [("Ni", 1, 40)]), "SiC": (3.21, [("Si", 1, 35), ("C", 1, 20)]),
    "SrTiO3": (5.12, [("Sr", 1, 80), ("Ti", 1, 70), ("O", 3, 45)]),
}


def published_srim(exe, runs, seeds, lines):
    lines += ["== Published SRIM damage calculations (Stoller et al. 2013 Table 4; Weber and Zhang 2019 Table 1)",
              "   AMCSET as TRIM with the papers' Ed, lattice binding 0, TRIM's Es; 1000 ions x 2 (light ions 20000);",
              "   T_dam: incident energy minus the ionization of ion and recoils (TRIM's energy partition)",
              f"   {'source':10s} {'case':22s} {'damage':6s} {'vac. paper':>10s} {'AMCSET':>17s} {'ratio':>6s}"
              f" {'T_dam paper':>11s} {'AMCSET':>9s} {'ratio':>6s}"]
    for src, ion, e_MeV, tkey, damage, vac, tdam in PUBLISHED_SRIM:
        rho, comp = PUBLISHED_TARGETS[tkey]
        v, td = [], []
        for k in range(2):
            sc = {"mode": "srim", "title": f"{ion} {e_MeV} MeV in {tkey}",
                  "ion": {"element": ion, "energy": f"{e_MeV} MeV"},
                  "layers": [{"name": tkey, "thickness": "1 um", "density": f"{rho} g/cm3",
                              "elements": [element(sym, ed, 0, None, f) for sym, f, ed in comp]}],
                  "ions": 20000 if ion == "H" else 1000, "seed": 101 + k, "threads": 0, "damage": damage,
                  "output": {"range_3d": False, "exit_lists": False}}
            s = run(exe, sc, runs / f"published_{ion}_{e_MeV:g}MeV_{tkey}_{damage}_s{k}", ("--auto-depth",))
            v.append(s["per_ion"]["vacancies"])
            ep = s["energy_percent"]
            td.append(e_MeV * 1000.0 * (100.0 - ep["ionization_ion"] - ep["ionization_recoils"]) / 100.0)
        (vm, ve), (tm, _) = mean_sem(v), mean_sem(td)
        record("published_srim", source=src, ion=ion, energy_MeV=e_MeV, target=tkey, damage=damage,
               paper_vacancies=vac, vacancies=vm, vacancies_sem=ve, paper_tdam_keV=tdam, tdam_keV=tm)
        tcol = f"{tdam:11.3f} {tm:9.3f} {tm / tdam:6.3f}" if tdam else f"{'-':>11s} {tm:9.3f} {'-':>6s}"
        lines.append(f"   {src:10s} {f'{ion} {e_MeV:g} MeV -> {tkey}':22s} {damage:6s} {vac:10.0f} {vm:9.1f} +- {ve:5.1f}"
                     f" {vm / vac:6.3f} {tcol}")
    lines.append("")


# A.F. Tasch, S.-H. Yang, S. Morris, D. Lim, J. Vac. Sci. Technol. B 12, 166 (1994),
# doi:10.1116/1.587177, Fig. 3 (values read off the figure): 15 keV As into (100) Si,
# 1e14 /cm2, rotation 0, SIMS. Peak about 6e19 /cm3 at about 12 nm; concentration at 100 and
# 200 nm: tilt 0: 3e17 and 3e16; tilt 10 deg: 1e16 and 2e15; amorphous Si: below 1e15 (the
# figure's floor) at both.
TASCH_AS = {"peak": 6e19, "0": (3e17, 3e16), "10": (1e16, 2e15), "amorphous": (1e15, 1e15)}
SI_CRYSTAL = {"structure": "diamond", "a": "5.431 A", "sublattices": ["Si"],
              "orientation": {"surface": [0, 0, 1], "reference": [1, 1, 0], "twist": "0 deg"},
              "temperature": "300 K", "vibrations": [{"element": "Si", "debye_temperature": "543 K"}]}


def as_si_channeling(exe, runs, seeds, lines):
    lines += ["== 15 keV As into (100) Si, channeling tails against tilt (Tasch et al. 1994, Fig. 3, read off the",
              "   figure). Concentration at 100 and 200 nm over the peak concentration. AMCSET: the crystal transport",
              "   (perfect crystal, 300 K, the default stopping division), 200000 ions; amorphous: trim2013.",
              f"   {'case':22s} {'peak depth nm':>13s} {'c(100)/peak':>16s} {'measured':>9s} {'c(200)/peak':>16s} {'measured':>9s}"]
    for label, tilt, crystal, oxide in (("amorphous Si", 0, False, False), ("crystal, tilt 0", 0, True, False),
                                        ("crystal, tilt 10 deg", 10, True, False),
                                        ("1 nm SiO2 + crystal, 0", 0, True, True),
                                        ("1 nm SiO2 + crystal, 10", 10, True, True)):
        layer = {"name": "Si", "thickness": "300 nm"}
        if crystal:
            layer["crystal"] = SI_CRYSTAL
        else:
            layer.update(density="2.329 g/cm3", elements=[element("Si")])
        layers = [layer]
        if oxide:  # the wafers' native oxide, about 1.0 nm (Tasch et al.)
            layers.insert(0, {"name": "SiO2", "thickness": "1 nm", "density": "2.32 g/cm3",
                              "elements": [element("Si", fraction=1), element("O", fraction=2)]})
        sc = {"mode": "srim", "title": f"As 15 keV into Si, {label}",
              "ion": {"element": "As", "energy": "15 keV", "angle": f"{tilt} deg"},
              "layers": layers, "ions": 200000, "seed": 101, "threads": 0, "damage": "quick",
              "physics": {"transport": "crystal" if crystal else "trim2013"},
              "output": {"range_3d": False, "exit_lists": False}}
        work = runs / f"tasch_as_{'crystal' if crystal else 'amorphous'}_{tilt}{'_oxide' if oxide else ''}"
        run(exe, sc, work)
        c, n, _ = ion_profile(work)
        sm = [sum(n[max(0, i - 1):i + 2]) / len(n[max(0, i - 1):i + 2]) for i in range(len(n))]
        peak = max(sm)
        pk = peak_depth(c, n) / 10.0

        def at(depth_nm):
            sel = [v for x, v in zip(c, n) if abs(x / 10.0 - depth_nm) <= 10.0]
            return sum(sel) / len(sel) / peak

        key = "amorphous" if not crystal else str(tilt)
        m100, m200 = (v / TASCH_AS["peak"] for v in TASCH_AS[key])
        fmt = (lambda v: f"< {v:7.1e}") if not crystal else (lambda v: f"{v:9.1e}")
        record("as_si_channeling", label=label, tilt_deg=tilt, crystal=crystal, oxide=oxide, peak_nm=pk,
               c100=at(100), c200=at(200), measured_c100=m100, measured_c200=m200,
               profile_depth_nm=[x / 10.0 for x in c], profile_over_peak=[v / peak for v in n])
        lines.append(f"   {label:22s} {pk:13.1f} {at(100):16.2e} {fmt(m100):>9s} {at(200):16.2e} {fmt(m200):>9s}")
    lines += ["   (AMCSET: mean over 90-110 and 190-210 nm, divided by the peak of the 3-bin smoothed profile;",
              "   the measured peak and depths are read off a log plot, about a factor 2 each)", ""]


def sb_si(exe, runs, seeds, lines):
    c = SB_SI
    vals, strag = [], []
    for k in range(seeds):
        sc = {"mode": "srim", "title": "Sb 50 keV in Si with native oxide",
              "ion": {"element": "Sb", "energy": f"{c['energy_keV']} keV", "angle": f"{c['angle_deg']} deg"},
              "layers": [{"name": "SiO2", "thickness": f"{c['oxide_nm']} nm", "density": "2.32 g/cm3",
                          "elements": [element("Si", 15, 2, 4.7, 1.0), element("O", 28, 3, 2, 2.0)]},
                         {"name": "Si", "thickness": "200 nm", "density": "2.329 g/cm3",
                          "elements": [element("Si", 15, 2, 4.7)]}],
              "ions": 2000, "seed": 101 + k, "threads": 0, "damage": "full"}
        s = run(exe, sc, runs / f"sb_si_s{k}")
        vals.append(s["range"]["mean_A"] / 10.0)
        strag.append(s["range"]["straggling_A"] / 10.0)
    m, e = mean_sem(vals)
    st, _ = mean_sem(strag)
    record("sb_si", measured_nm=c["measured_nm"], measured_err_nm=c["measured_err_nm"],
           measured_straggle_nm=c["straggle_nm"], paper_srim_nm=c["srim_nm"], amcset_nm=m, amcset_sem=e,
           amcset_straggling_nm=st)
    lines += ["== 50 keV Sb into Si (7 deg, 2.9 nm native oxide): mean depth (Titze et al. 2021, SIMS)",
              f"   measured {c['measured_nm']} +- {c['measured_err_nm']} nm (straggle {c['straggle_nm']} nm); "
              f"their SRIM-2013 run {c['srim_nm']} nm",
              f"   AMCSET (full cascades, 2000 ions x {seeds}): {m:.2f} +- {e:.2f} nm (straggling {st:.2f} nm); "
              f"AMCSET / measured {m / c['measured_nm']:.3f}, AMCSET / their SRIM {m / c['srim_nm']:.3f}", ""]


def isotropic(rng):
    cth = 2.0 * rng.random() - 1.0
    phi = 2.0 * math.pi * rng.random()
    s = math.sqrt(max(0.0, 1.0 - cth * cth))
    return [cth, s * math.cos(phi), s * math.sin(phi)]


def recoil_yield(exe, runs, sym, rho, e_keV, seeds, pkas=500, model="srim"):
    """Ionization share of PKAs of energy e_keV deep in a pure target (full cascades),
    with the electronic stopping model `model` (srim: the bundled SRIM-2013 tables;
    lss: Lindhard-Scharff)."""
    vals = []
    for k in range(seeds):
        rng = random.Random(f"{sym}-{e_keV}-{k}")
        particles = [{"element": sym, "energy": f"{e_keV * 1000.0:.6g} eV", "x": "5000 A",
                      "direction": isotropic(rng)} for _ in range(pkas)]
        sc = {"mode": "srim", "title": f"{sym} PKA {e_keV} keV",
              "ion": {"element": sym, "energy": f"{e_keV} keV"},
              "layers": [{"name": sym, "thickness": "1 um", "density": f"{rho} g/cm3", "elements": [element(sym)]}],
              "source": {"type": "recoils", "particles": particles}, "seed": 101 + k, "threads": 0,
              "damage": "full", "output": {"range_3d": False, "exit_lists": False}}
        tag = "" if model == "srim" else f"_{model}"
        s = run(exe, sc, runs / f"yield_{sym}_{e_keV:g}keV{tag}_s{k}", ("--stopping", model))
        ep = s["energy_percent"]
        vals.append((ep["ionization_ion"] + ep["ionization_recoils"]) / 100.0)
    return mean_sem(vals)


def quenching(exe, runs, seeds, lines):
    lines += ["== Ionization yield of nuclear recoils (ionization / recoil energy); AMCSET: PKAs inside a",
              "   pure target, full cascades, 500 PKAs x seeds, TRIM's default Ed/Eb/Es; Lindhard: k g/(1 + k g)",
              "   Ge in Ge (Bonhomme et al. 2022; Lindhard k = 0.162, their fit)",
              f"   {'E (keV)':>8s} {'measured':>16s} {'AMCSET':>16s} {'AMCSET, lss':>16s} {'Lindhard':>9s}"]
    for e, q, dq in GE_YIELD:
        m, s = recoil_yield(exe, runs, "Ge", 5.323, e, seeds)
        ml, sl = recoil_yield(exe, runs, "Ge", 5.323, e, seeds, model="lss")
        record("quenching", target="Ge", source="Bonhomme 2022", energy_keV=e, measured=q, measured_err=dq,
               amcset=m, amcset_lss=ml, lindhard=lindhard_yield(e, 32, 0.162))
        lines.append(f"   {e:8.2f} {q:8.3f} +- {dq:.3f} {m:8.3f} +- {s:.3f} {ml:8.3f} +- {sl:.3f} "
                     f"{lindhard_yield(e, 32, 0.162):9.3f}")
    lines += ["   Si in Si (Chavarria et al. 2016: E_e/E_r with the uncertainty of E_r; Izraelevitch et al. 2017;",
              "   Lindhard k = 0.15, the value Chavarria et al. compare with)",
              f"   {'E (keV)':>8s} {'measured':>16s} {'AMCSET':>16s} {'AMCSET, lss':>16s} {'Lindhard':>9s}"]
    for ee, er, der in SI_CHAVARRIA:
        q = ee / er
        dq = q * der / er
        m, s = recoil_yield(exe, runs, "Si", 2.329, er, seeds)
        ml, sl = recoil_yield(exe, runs, "Si", 2.329, er, seeds, model="lss")
        record("quenching", target="Si", source="Chavarria 2016", energy_keV=er, measured=q, measured_err=dq,
               amcset=m, amcset_lss=ml, lindhard=lindhard_yield(er, 14, 0.15))
        lines.append(f"   {er:8.2f} {q:8.3f} +- {dq:.3f} {m:8.3f} +- {s:.3f} {ml:8.3f} +- {sl:.3f} "
                     f"{lindhard_yield(er, 14, 0.15):9.3f}  (C)")
    for er, q, dq in SI_IZRAELEVITCH:
        m, s = recoil_yield(exe, runs, "Si", 2.329, er, seeds)
        ml, sl = recoil_yield(exe, runs, "Si", 2.329, er, seeds, model="lss")
        record("quenching", target="Si", source="Izraelevitch 2017", energy_keV=er, measured=q, measured_err=dq,
               amcset=m, amcset_lss=ml, lindhard=lindhard_yield(er, 14, 0.15))
        lines.append(f"   {er:8.2f} {q:8.3f} +- {dq:.3f} {m:8.3f} +- {s:.3f} {ml:8.3f} +- {sl:.3f} "
                     f"{lindhard_yield(er, 14, 0.15):9.3f}  (I)")
    lines.append("")


CASES = {"ni_ranges": ni_ranges, "ni_peaks": ni_peaks, "ranges": ranges, "sputtering": sputtering,
         "published_srim": published_srim, "as_si_channeling": as_si_channeling, "sb_si": sb_si,
         "quenching": quenching}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--only", default=",".join(CASES))
    ap.add_argument("--seeds", type=int, default=4)
    a = ap.parse_args()
    runs = a.out / "runs"
    lines = [f"AMCSET against measurements (tools/validation/literature_cases.py, {a.seeds} seeds)", ""]
    for name in a.only.split(","):
        CASES[name](a.exe, runs, a.seeds, lines)
    a.out.mkdir(parents=True, exist_ok=True)
    (a.out / "results.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    old = {}
    if (a.out / "results.json").exists():  # a partial rerun (--only) keeps the other cases
        old = json.loads((a.out / "results.json").read_text(encoding="utf-8"))
    old.update(RESULTS)
    (a.out / "results.json").write_text(json.dumps(old, indent=1), encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

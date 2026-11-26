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

"""The figures of the validation review (docs/validation/), drawn from the output of
the validation scripts. Unlike the other tools this one needs matplotlib
(`pip install matplotlib`); it only draws, it computes no evidence.

    python tools/validation/make_figures.py [--exe build/Release/amcset.exe] [--only NAME,...]

Inputs (defaults relative to the repository root; each figure is skipped with a
message when its input is missing):
  build/validation/literature/   literature_cases.py: results.json, ranges.csv, sputtering.csv
  build/validation/stopping/     stopping_reference.py: iaea_ratios.csv
  build/validation/table_scan/   srim_table_scan.py: ratios.csv
  build/validation/sources/data/ NIST pstar_*.txt, astar_*.txt
  build/trim_validation/         trim3/analysis.json, trim4/analysis.json, trim4_nickel/fix_effect.txt
  build/damage_models/           eb0/damage_functions.csv, mitsi2023/comparison.txt
  build/crystal_calibration/     results.txt
and `amcset stopping` (--exe, with its DLLs on PATH) for stopping curves.
Output: docs/validation/figures/<name>.svg
"""

import argparse
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path

import logging

import matplotlib

matplotlib.use("svg")
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "trim_compare"))

# --- style: the reference palette of the review ---------------------------------------
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = C
RAMP = {250: "#86b6ef", 450: "#2a78d6", 650: "#104281"}  # ordinal blue steps
SURF, INK, INK2, MUTED, GRID, AXIS, BAND = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#f0efec"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 9.5, "axes.titleweight": "semibold",
    "axes.titlelocation": "left", "axes.prop_cycle": matplotlib.cycler(color=C),
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "xtick.major.size": 3, "ytick.major.size": 3, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 8.5,
    "legend.frameon": False, "legend.fontsize": 8, "legend.labelcolor": INK2,
    "lines.linewidth": 1.5, "lines.markersize": 5.5, "lines.markeredgecolor": SURF, "lines.markeredgewidth": 0.8,
    "errorbar.capsize": 0, "svg.fonttype": "none", "svg.hashsalt": "amcset", "figure.dpi": 100,
})

OUT = ROOT / "docs" / "validation" / "figures"
PNG = None  # --png DIR: also write PNG previews there
FIGURES = {}


def figure(name):
    def wrap(f):
        FIGURES[name] = f
        return f
    return wrap


def log_label(v, _pos=None):
    if v <= 0:
        return ""
    if 1e-4 <= v < 1e5:
        return f"{v:.10g}"
    e = int(math.floor(math.log10(v) + 1e-9))
    m = round(v / 10 ** e, 6)
    return f"$10^{{{e}}}$" if m == 1 else f"${m:g}\\times10^{{{e}}}$"


def tidy_log(ax):
    """Log axes: plain numbers, no minor labels, 1-2-5 ticks on short spans."""
    from matplotlib.ticker import FixedFormatter, FixedLocator, FuncFormatter, LogLocator, NullFormatter
    for axis, lim, scale in ((ax.xaxis, ax.get_xlim(), ax.get_xscale()), (ax.yaxis, ax.get_ylim(), ax.get_yscale())):
        if scale != "log":
            continue
        axis.set_minor_formatter(NullFormatter())
        if isinstance(axis.get_major_locator(), FixedLocator) or isinstance(axis.get_major_formatter(), FixedFormatter):
            continue
        lo, hi = sorted(lim)
        if lo > 0 and math.log10(hi / lo) < 2.5:
            axis.set_major_locator(LogLocator(base=10, subs=(1, 2, 5)))
        axis.set_major_formatter(FuncFormatter(log_label))


def save(fig, name, caption=None):
    OUT.mkdir(parents=True, exist_ok=True)
    for ax in fig.axes:
        tidy_log(ax)
    if caption:  # below everything the axes draw (tick labels, axis labels)
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        r = FigureCanvasAgg(fig).get_renderer()
        fig.set_canvas(fig.canvas)
        bottom = min(ax.get_tightbbox(r).y0 for ax in fig.axes) / fig.bbox.height
        fig.text(0.0, bottom - 0.015, caption, ha="left", va="top", fontsize=7.5, color=MUTED, wrap=True)
    fig.savefig(OUT / f"{name}.svg", bbox_inches="tight", metadata={"Date": None, "Creator": "AMCSET make_figures.py"})
    if PNG:
        PNG.mkdir(parents=True, exist_ok=True)
        fig.savefig(PNG / f"{name}.png", bbox_inches="tight", dpi=110)
    plt.close(fig)
    print("wrote", name)


def ratio_axes(ax, lo=None, hi=None, band=0.1, log=False):
    """Reference band (1 +- band) and the line at 1 on a ratio axis."""
    ax.axhspan(1 - band, 1 + band, color=BAND, zorder=0, lw=0)
    ax.axhline(1.0, color=AXIS, lw=0.9, zorder=1)
    if log:
        ax.set_yscale("log")
    if lo is not None:
        ax.set_ylim(lo, hi)


def vratio_axes(ax, band=0.1):
    ax.axvspan(1 - band, 1 + band, color=BAND, zorder=0, lw=0)
    ax.axvline(1.0, color=AXIS, lw=0.9, zorder=1)


def meas_style(**kw):
    """Measured data: ink-coloured markers."""
    d = dict(color=INK, marker="o", ls="none", ms=4.5, mec=SURF, mew=0.6, zorder=5)
    d.update(kw)
    return d


def read_csv(path):
    with open(path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


class Inputs:
    def __init__(self, a):
        self.a = a
        self._cache = {}

    def need(self, path):
        if not Path(path).exists():
            raise FileNotFoundError(path)
        return Path(path)

    def lit(self):
        if "lit" not in self._cache:
            self._cache["lit"] = json.loads(self.need(self.a.literature / "results.json").read_text(encoding="utf-8"))
        return self._cache["lit"]

    def ranges(self):
        return read_csv(self.need(self.a.literature / "ranges.csv"))

    def sputter(self):
        return read_csv(self.need(self.a.literature / "sputtering.csv"))

    def iaea(self):
        if "iaea" not in self._cache:
            self._cache["iaea"] = read_csv(self.need(self.a.stopping / "iaea_ratios.csv"))
        return self._cache["iaea"]

    def scan(self):
        return read_csv(self.need(self.a.scan / "ratios.csv"))

    def trim(self, name):
        return json.loads(self.need(self.a.trim / name / "analysis.json").read_text(encoding="utf-8"))

    def table(self, ion, target, lo, hi, model="srim", fix=False, density=None):
        import stopping_reference as sr
        if not self.a.exe:
            raise FileNotFoundError("--exe (amcset) for stopping curves")
        key = (ion, target, lo, hi, model, fix, density)
        if key not in self._cache:
            self._cache[key] = sr.table(self.a.exe, ion, target, lo, hi, model, False, fix, density)
        return self._cache[key]

    def mass(self, ion):
        import stopping_reference as sr
        return sr.default_mass(self.a.exe, ion)


# =========================================================================================
# 0. Overview
# =========================================================================================

def pct(v, p):
    v = sorted(v)
    k = (len(v) - 1) * p
    i = int(math.floor(k))
    j = min(i + 1, len(v) - 1)
    return v[i] + (v[j] - v[i]) * (k - i)


@figure("overview")
def fig_overview(inp):
    rows = []  # (label, group, values)
    tr = inp.trim("trim3")
    import doc_tables
    vals = []
    for q, _ in doc_tables.QUANTITIES:
        for u in tr.values():
            v = u["quantities"].get(q)
            if v and v.get("diff_percent") is not None and v["trim"]:
                vals.append(v["amcset"] / v["trim"])
    rows.append(("TRIM.exe: 1189 quantities, 75 inputs", "Against TRIM.exe", vals))
    pub = inp.lit()["published_srim"]
    rows.append(("Published SRIM runs (vacancies)", "Against TRIM.exe", [r["vacancies"] / r["paper_vacancies"] for r in pub]))
    ia = inp.iaea()
    rows.append(("IAEA stopping, 38 268 measurements", "Stopping", [float(r["srim"]) / float(r["measured"]) for r in ia]))
    ni = inp.lit()["ni_peaks"]
    rows.append(("Ni peak depths: SRIM-2013 tables", "Ni ions", [r["srim_um"] / r["measured_um"] for r in ni]))
    rows.append(("Ni peak depths: nickel_ions", "Ni ions", [r["fix_um"] / r["measured_um"] for r in ni]))
    rg = inp.ranges()
    for label, groups in (("Ranges in Si, SiO2 (193)", ("behar_si", "behar_sio2", "behar_light", "oosterhoff")),
                          ("Ranges in C, B, SiC (157)", ("behar_c", "behar_b", "behar_sic", "zhang_sic"))):
        rows.append((label, "Ranges", [float(r["ratio"]) for r in rg if r["group"] in groups]))
    sp = [r for r in inp.sputter() if not r["flag"] and float(r["energy_keV"]) >= 0.5]
    rows.append(("Sputtering yields, E >= 0.5 keV (98)", "Sputtering", [float(r["amcset"]) / float(r["measured"]) for r in sp]))
    q = inp.lit()["quenching"]
    rows.append(("Ionization yield of keV recoils (19)", "Ionization", [r["amcset"] / r["measured"] for r in q]))

    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    vratio_axes(ax)
    for i, (label, _, v) in enumerate(rows):
        y = len(rows) - 1 - i
        lo, q1, med, q3, hi = pct(v, 0.05), pct(v, 0.25), statistics.median(v), pct(v, 0.75), pct(v, 0.95)
        ax.plot([lo, hi], [y, y], color=AXIS, lw=1.2, solid_capstyle="round", zorder=2)
        ax.plot([q1, q3], [y, y], color=BLUE, lw=4.5, solid_capstyle="round", zorder=3)
        ax.plot(med, y, "o", color=INK, ms=5, zorder=4)
        ax.text(hi * 1.03, y, f"{med:.3f}", va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r[0] for r in rows][::-1])
    ax.set_xscale("log")
    ax.set_xlim(0.35, 3.2)
    ax.set_xticks([0.4, 0.5, 0.7, 1, 1.4, 2, 3])
    ax.set_xticklabels(["0.4", "0.5", "0.7", "1", "1.4", "2", "3"])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("AMCSET / reference")
    ax.set_title("Every comparison at a glance: median (dot), middle half (bar), 5-95 % (line)")
    save(fig, "overview", "Shaded: within 10 %. Each row is one study of the review; the number is the median ratio.")


# =========================================================================================
# 1. TRIM.exe
# =========================================================================================

@figure("trim_z")
def fig_trim_z(inp):
    import doc_tables
    tr = inp.trim("trim3")
    zs = [z for q, _ in doc_tables.QUANTITIES for z, _ in doc_tables.z_values(tr, q)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.6), gridspec_kw={"width_ratios": [1, 1.25]})
    bins = [x * 0.5 - 4.0 for x in range(17)]
    a.hist(zs, bins=bins, color=BLUE, edgecolor=SURF, linewidth=1.5)
    xs = [x / 20.0 - 4.0 for x in range(161)]
    a.plot(xs, [len(zs) * 0.5 * math.exp(-x * x / 2) / math.sqrt(2 * math.pi) for x in xs], color=INK, lw=1.2)
    a.set_xlabel("z  (AMCSET − TRIM.exe, in standard errors)")
    a.set_ylabel("comparisons")
    a.set_title(f"All {len(zs)} comparisons against a normal distribution")
    a.text(-3.9, a.get_ylim()[1] * 0.9, "line: normal distribution\nN(0, 1) for the same count", color=INK2, fontsize=7.5,
           va="top")
    labels = []
    for i, (q, label) in enumerate(doc_tables.QUANTITIES):
        z = [v for v, _ in doc_tables.z_values(tr, q)]
        y = len(doc_tables.QUANTITIES) - 1 - i
        b.plot(z, [y] * len(z), "o", color=BLUE, ms=3.6, alpha=0.55, mew=0)
        labels.append(label)
    b.axvspan(-2, 2, color=BAND, zorder=0, lw=0)
    b.axvline(0, color=AXIS, lw=0.9)
    for x in (-3, 3):
        b.axvline(x, color=MUTED, lw=0.6)
    b.set_yticks(range(len(labels)))
    b.set_yticklabels(labels[::-1], fontsize=7.5)
    b.grid(axis="y", visible=False)
    b.set_xlabel("z")
    b.set_title("Per quantity (shaded |z| < 2, lines |z| = 3)")
    save(fig, "trim_z", "161 TRIM.exe runs of 75 inputs against 704 AMCSET runs (build/trim_validation/trim3/analysis.json). "
         "The one point beyond 3 is replacements in Xe gas.")


@figure("trim_profiles")
def fig_trim_profiles(inp):
    tr = inp.trim("trim3")
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    bins = [0.4 + 0.1 * i for i in range(20)]
    for key, color, label in (("vacancies", BLUE, "vacancy profiles"), ("ion ranges", ORANGE, "ion range profiles")):
        v = [u["profiles"][key]["chi2_per_ndf"] for u in tr.values() if (u.get("profiles") or {}).get(key)]
        ax.hist(v, bins=bins, histtype="step", color=color, lw=1.8, label=f"{label} ({len(v)}, median {statistics.median(v):.2f})")
    ax.axvline(1.0, color=AXIS, lw=0.9)
    ax.set_xlabel("χ² per degree of freedom (AMCSET's profile against TRIM.exe's, bin by bin)")
    ax.set_ylabel("inputs")
    ax.legend(loc="upper right")
    ax.set_title("Depth profiles against TRIM.exe")
    save(fig, "trim_profiles")


@figure("trim_nickel_batch")
def fig_trim_nickel(inp):
    path = inp.need(inp.a.trim / "trim4_nickel" / "fix_effect.txt")
    rows = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"(\S+_kp|\S+_fc)\s+(.+?)\s{2,}(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)$", line)
        if m:
            case, q, t, a, z, f, r = m.groups()
            rows.setdefault(case, {})[q] = (float(t), float(a), float(z), float(f))
    quants = [("mean projected range (A)", "Mean projected range"), ("range straggling (A)", "Range straggling"),
              ("vacancies / ion", "Vacancies per ion"), ("ionization by ion (%)", "Ionization by the ion")]
    order = sorted(rows, key=lambda c: (c.split("_")[1], c.split("_")[0], int(re.search(r"_(\d+)keV", c).group(1))))
    fig, axes = plt.subplots(1, 4, figsize=(10.5, 3.8), sharey=True)
    for ax, (q, title) in zip(axes, quants):
        vratio_axes(ax, band=0.02)
        for i, c in enumerate(order):
            if q not in rows[c]:
                continue
            t, a, z, f = rows[c][q]
            y = len(order) - 1 - i
            ax.plot(a / t, y, "o", color=BLUE)
            ax.plot(f / t, y, "D", color=ORANGE)
        ax.set_title(title, fontsize=8.5)
        ax.set_xlabel("AMCSET / TRIM.exe")
        ax.grid(axis="y", visible=False)
    def case_label(c):
        ion, target, e, mode = re.match(r"(\w+?)_(\w+?)_(\d+)keV_(\w+)", c).groups()
        return f"{ion.capitalize()} → {target.capitalize()}, {int(e) / 1000:g} MeV ({mode})"
    axes[0].set_yticks(range(len(order)))
    axes[0].set_yticklabels([case_label(c) for c in order][::-1], fontsize=7.5)
    fig.legend([Line2D([], [], marker="o", color=BLUE, ls="none"), Line2D([], [], marker="D", color=ORANGE, ls="none")],
               ["AMCSET as TRIM", "with nickel_ions"], loc="upper right", bbox_to_anchor=(1.0, 1.06), ncol=2)
    fig.suptitle("Ni, Co and Cu ions: AMCSET against TRIM.exe, without and with the Ni correction (batch 4)",
                 x=0.0, y=1.03, ha="left", fontsize=9.5, fontweight="semibold")
    save(fig, "trim_nickel_batch", "Shaded: within 2 %. kp: quick calculation, fc: full cascades. One TRIM run per case, 8 AMCSET seeds.")


@figure("published_srim")
def fig_published(inp):
    pub = inp.lit()["published_srim"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    vratio_axes(ax, band=0.05)
    for i, r in enumerate(pub):
        y = len(pub) - 1 - i
        color = BLUE if r["source"].startswith("Stoller") else ORANGE
        ax.plot(r["vacancies"] / r["paper_vacancies"], y, "o", color=color, ms=6)
        if r["paper_tdam_keV"]:
            ax.plot(r["tdam_keV"] / r["paper_tdam_keV"], y, "s", color=color, mfc=SURF, mec=color, mew=1.2, ms=5)
    ax.set_yticks(range(len(pub)))
    ax.set_yticklabels([f"{r['ion']} {r['energy_MeV']:g} MeV → {r['target']} ({r['damage']})" for r in pub][::-1], fontsize=7.5)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("AMCSET / the paper's SRIM")
    ax.legend([Line2D([], [], marker="o", color=BLUE, ls="none"), Line2D([], [], marker="o", color=ORANGE, ls="none"),
               Line2D([], [], marker="s", color=INK2, mfc=SURF, ls="none")],
              ["Stoller et al. 2013", "Weber and Zhang 2019", "damage energy (Stoller)"], loc="lower left", fontsize=7.5)
    ax.set_title("Published SRIM damage calculations: vacancies per ion")
    save(fig, "published_srim", "Shaded: within 5 %. The proton's damage energy is defined differently in the paper (its Eq. 5b).")


# =========================================================================================
# 2. Stopping
# =========================================================================================

@figure("iaea_histogram")
def fig_iaea_hist(inp):
    ia = inp.iaea()
    fig, ax = plt.subplots(figsize=(6.2, 3.3))
    bins = [10 ** (-0.4 + 0.02 * i) for i in range(41)]
    for key, color, label in (("srim", BLUE, "SRIM-2013 tables (default)"), ("builtin", ORANGE, "builtin (no tables)")):
        v = [float(r[key]) / float(r["measured"]) for r in ia]
        ax.hist(v, bins=bins, histtype="step", color=color, lw=1.8,
                label=f"{label}: median {statistics.median(v):.3f}")
    ax.set_xscale("log")
    ax.axvspan(0.9, 1.1, color=BAND, zorder=0, lw=0)
    ax.set_xticks([0.4, 0.5, 0.7, 1, 1.4, 2, 2.5])
    ax.set_xticklabels(["0.4", "0.5", "0.7", "1", "1.4", "2", "2.5"])
    ax.set_xlabel("AMCSET electronic stopping / measured")
    ax.set_ylabel("measurements")
    ax.legend(loc="upper left")
    ax.set_title(f"{len(ia)} measurements in elemental targets (IAEA database)")
    save(fig, "iaea_histogram", "Shaded: within 10 %. Calculated rows of the database (SRIM-2013, MSTAR, DPASS, ...) left out.")


def running_median(xs, ys, nbins=30):
    pts = sorted(zip(xs, ys))
    lo, hi = math.log10(pts[0][0]), math.log10(pts[-1][0])
    out = []
    for i in range(nbins):
        a, b = lo + (hi - lo) * i / nbins, lo + (hi - lo) * (i + 1) / nbins
        sel = [y for x, y in pts if a <= math.log10(x) < b]
        if len(sel) >= 15:
            out.append((10 ** ((a + b) / 2), statistics.median(sel), pct(sel, 0.16), pct(sel, 0.84)))
    return out


@figure("iaea_energy")
def fig_iaea_energy(inp):
    ia = inp.iaea()
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, key, title in zip(axes, ("srim", "builtin"), ("SRIM-2013 tables", "builtin stopping (no tables)")):
        x = [float(r["keV_per_u"]) / 1000.0 for r in ia]
        y = [float(r[key]) / float(r["measured"]) for r in ia]
        ratio_axes(ax, 0.3, 3.0, log=True)
        ax.plot(x, y, ".", color=BLUE, ms=1.6, alpha=0.18, mew=0, rasterized=True)
        rm = running_median(x, y)
        ax.fill_between([r[0] for r in rm], [r[2] for r in rm], [r[3] for r in rm], color=ORANGE, alpha=0.18, lw=0)
        ax.plot([r[0] for r in rm], [r[1] for r in rm], color=ORANGE, lw=1.8)
        ax.set_xscale("log")
        ax.set_xlim(1e-4, 2e3)
        ax.set_xlabel("energy per mass (MeV/u)")
        ax.set_title(title)
        ax.set_yticks([0.3, 0.5, 0.7, 1, 1.4, 2, 3])
        ax.set_yticklabels(["0.3", "0.5", "0.7", "1", "1.4", "2", "3"])
    axes[0].set_ylabel("AMCSET / measured")
    axes[1].legend([Line2D([], [], marker=".", color=BLUE, ls="none"), Line2D([], [], color=ORANGE)],
                   ["one measurement", "running median, 16-84 %"], loc="upper right")
    save(fig, "iaea_energy", "IAEA electronic stopping power database, release 2026-09; every ion in every elemental target.")


@figure("iaea_per_ion")
def fig_iaea_per_ion(inp):
    import stopping_reference as sr
    ia = inp.iaea()
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 5.6), sharex=True, sharey=True)
    for ax, (label, a, b) in zip(axes.flat, sr.BANDS):
        by = {}
        for r in ia:
            e = float(r["keV_per_u"])
            if a <= e < b:
                by.setdefault(r["ion"], []).append(float(r["srim"]) / float(r["measured"]))
        ratio_axes(ax, 0.6, 1.5)
        zs = [(sr.Z[i], statistics.median(v), len(v), i) for i, v in by.items() if len(v) >= 10]
        ax.plot([z[0] for z in zs], [z[1] for z in zs], "o", color=BLUE, ms=4.5)
        for z, m, n, ion in zs:
            if ion in ("Ni", "Au") or abs(m - 1) > 0.12:
                ax.annotate(ion, (z, m), xytext=(3, 3), textcoords="offset points", fontsize=7, color=INK2)
        ax.set_title(f"{label} MeV/u" if "MeV" not in label else label)
    for ax in axes[1]:
        ax.set_xlabel("ion Z1")
    for ax in axes[:, 0]:
        ax.set_ylabel("median AMCSET / measured")
    save(fig, "iaea_per_ion", "Median per ion (at least 10 measurements), bundled SRIM-2013 tables; labelled: Ni, Au and medians off by more than 12 %.")


@figure("ni_stopping_in_fe")
def fig_ni_stopping(inp):
    fig, ax = plt.subplots(figsize=(6.6, 3.8))
    for ion in ("Mn", "Fe", "Co", "Cu", "Zn"):
        m = inp.mass(ion)
        t = inp.table(ion, "Fe", 5.0 * m, 10000.0 * m)
        xs = [r[0] / m / 1000.0 for r in t]
        ax.plot(xs, [r[1] for r in t], color=MUTED, lw=0.9, label="Mn, Fe, Co, Cu, Zn" if ion == "Mn" else None)
    m = inp.mass("Ni")
    t = inp.table("Ni", "Fe", 5.0 * m, 10000.0 * m)
    f = inp.table("Ni", "Fe", 5.0 * m, 10000.0 * m, fix=True)
    ax.plot([r[0] / m / 1000.0 for r in t], [r[1] for r in t], color=ORANGE, lw=2.0, label="Ni, SRIM-2013 tables")
    ax.plot([r[0] / m / 1000.0 for r in f], [r[1] for r in f], color=BLUE, lw=2.0, label="Ni with nickel_ions")
    ax.axvspan(0.02, 0.3, color=BAND, zorder=0, lw=0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, 10)
    ax.set_xlabel("energy per mass (MeV/u)")
    ax.set_ylabel("electronic stopping in Fe (MeV cm²/mg)")
    ax.legend(loc="upper left")
    ax.set_title("Ni ions and their neighbours Mn-Zn in iron (bundled tables)")
    save(fig, "ni_stopping_in_fe", "Shaded: 0.02-0.3 MeV/u, where SRIM-2013's Ni stopping dips below its neighbours'. Grey: Mn, Fe, Co, Cu, Zn.")


@figure("table_scan")
def fig_table_scan(inp):
    rows = [r for r in inp.scan() if r["target"] == "Fe"]
    es = sorted({float(r["MeV_per_u"]) for r in rows})
    e0 = min(es, key=lambda e: abs(math.log(e / 0.1)))
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.5, 3.5), gridspec_kw={"width_ratios": [1.3, 1]})
    sel = [r for r in rows if float(r["MeV_per_u"]) == e0]
    ratio_axes(a, 0.5, 1.6, band=0.15)
    for r in sel:
        z, v = int(r["z1"]), fnum(r["ratio_z1pm1"])
        if v is None:
            continue
        hit = r["ion"] in ("Ni", "Au")
        a.plot(z, v, "o", color=ORANGE if hit else BLUE, ms=5.5 if hit else 4)
        if hit:
            a.annotate(r["ion"], (z, v), xytext=(4, -10), textcoords="offset points", fontsize=8, color=INK)
    a.set_xlabel("ion Z1")
    a.set_ylabel("S_e / geometric mean of Z1 ± 1")
    a.set_title(f"Every ion in Fe at {e0:.3g} MeV/u")
    ratio_axes(b, 0.5, 1.2, band=0.15)
    for ion, color in (("Ni", ORANGE), ("Au", VIOLET)):
        s = sorted((float(r["MeV_per_u"]), fnum(r["ratio_z1pm1"])) for r in rows if r["ion"] == ion)
        b.plot([x for x, _ in s], [y for _, y in s], color=color, lw=1.8, label=ion)
        s2 = sorted((float(r["MeV_per_u"]), fnum(r["ratio_z1pm2"])) for r in rows if r["ion"] == ion)
        b.plot([x for x, _ in s2], [y for _, y in s2], color=color, lw=1.0, ls=(0, (1, 1.5)))
    b.set_xscale("log")
    b.set_xlabel("energy per mass (MeV/u)")
    b.legend(loc="lower right", title="solid: Z1 ± 1, dotted: Z1 ± 2", title_fontsize=7)
    b.set_title("Ni and Au against their neighbours in Fe")
    save(fig, "table_scan", "Scan of the bundled SRIM-2013 tables (srim_table_scan.py). Shaded: 1 ± 0.15. The ratios are the same in all nine targets scanned.")


def ni_dataset(ref):
    if ref.startswith("Zha2002"):
        return 0
    if ref.startswith("Ay198"):
        return 1
    return 2


@figure("iaea_nickel")
def fig_iaea_ni(inp):
    ia = [r for r in inp.iaea() if r["ion"] == "Ni"]
    names = ["Zhang 2002 (ToF-ERDA)", "Anthony and Lanford 1981-82", "other data sets"]
    marks = ["o", "s", "^"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, key, title in zip(axes, ("srim", "srim_nickel_fix"), ("SRIM-2013 tables", "with nickel_ions")):
        ratio_axes(ax, 0.5, 2.0, log=True)
        ax.axvspan(0.02, 0.3, color="#eef4fc", zorder=0, lw=0)
        for k in range(3):
            sel = [r for r in ia if ni_dataset(r["ref"]) == k]
            ax.plot([float(r["keV_per_u"]) / 1000 for r in sel], [float(r[key]) / float(r["measured"]) for r in sel],
                    marks[k], color=C[k], ms=4.5, label=names[k])
        ax.set_xscale("log")
        ax.set_xlabel("energy per mass (MeV/u)")
        ax.set_title(title)
        ax.set_yticks([0.5, 0.7, 1, 1.4, 2])
        ax.set_yticklabels(["0.5", "0.7", "1", "1.4", "2"])
    axes[0].set_ylabel("AMCSET / measured, Ni ions")
    axes[0].legend(loc="upper left", fontsize=7.5)
    save(fig, "iaea_nickel", "Every measured Ni-ion stopping in an elemental target (C, Al, Ti, Ni, Cu, Ag, Ta, Au, Ar). Blue band: SRIM's Ni dip; grey: within 10 %.")


@figure("zhang2002_neighbours")
def fig_zhang2002(inp):
    ia = [r for r in inp.iaea() if r["ref"].startswith("Zha2002") and r["ion"] in ("Fe", "Co", "Ni", "Cu")]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.4), sharey=True)
    for ax, target in zip(axes, ("C", "Al")):
        ratio_axes(ax, 0.7, 1.25)
        nb = [r for r in ia if r["target"] == target and r["ion"] != "Ni"]
        ni = [r for r in ia if r["target"] == target and r["ion"] == "Ni"]
        ax.plot([float(r["keV_per_u"]) for r in nb], [float(r["srim"]) / float(r["measured"]) for r in nb], "o",
                color=MUTED, ms=4, label="Fe, Co, Cu (SRIM tables)")
        ax.plot([float(r["keV_per_u"]) for r in ni], [float(r["srim"]) / float(r["measured"]) for r in ni], "o",
                color=ORANGE, ms=5.5, label="Ni, SRIM tables")
        ax.plot([float(r["keV_per_u"]) for r in ni], [float(r["srim_nickel_fix"]) / float(r["measured"]) for r in ni],
                "D", color=BLUE, ms=5, label="Ni with nickel_ions")
        ax.set_title(f"in {'carbon' if target == 'C' else 'aluminium'}")
        ax.set_xlabel("energy per mass (keV/u)")
    axes[0].set_ylabel("AMCSET / measured")
    axes[1].legend(loc="lower right", fontsize=7.5)
    save(fig, "zhang2002_neighbours", "Y. Zhang, Nucl. Instrum. Methods B 196, 1 (2002): Fe, Co, Ni and Cu measured in the same apparatus (from the IAEA database).")


def nist_rows(path):
    import stopping_reference as sr
    return sr.read_nist(path)


NIST_MATS = [("silicon", "Si", None), ("iron", "Fe", None), ("copper", "Cu", None), ("tungsten", "W", None),
             ("gold", "Au", None), ("water", "H:2,O:1", "1.0")]


def nist_figure(inp, prog, ion, mass, name, title):
    fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.4), sharex=True, sharey=True)
    for ax, (mat, target, density) in zip(axes.flat, NIST_MATS):
        ref = nist_rows(inp.need(inp.a.nist / f"{prog}_{mat}.txt"))
        lo, hi = ref[0][0] * 1000 / 1.2, ref[-1][0] * 1000 * 1.2
        ratio_axes(ax, 0.6, 1.6, band=0.05)
        for model, color, label in (("srim", BLUE, "SRIM-2013 tables"), ("builtin", ORANGE, "builtin"),
                                    ("bethe-bloch", AQUA, "Bethe-Bloch")):
            t = inp.table(ion, target, lo, hi, model, density=density)
            import stopping_reference as sr
            pts = [(r[0] / mass, sr.interpolate(t, r[0] * 1000.0, 1) * 1000.0 / r[1]) for r in ref]
            if model == "bethe-bloch":
                pts = [p for p in pts if p[0] >= 1.0]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=1.6, label=label)
        ax.set_xscale("log")
        ax.set_title(mat)
    for ax in axes[1]:
        ax.set_xlabel("energy per mass (MeV/u)")
    for ax in axes[:, 0]:
        ax.set_ylabel("AMCSET / NIST, electronic")
    axes[0, 0].legend(loc="upper right", fontsize=7)
    fig.suptitle(title, x=0.0, y=1.03, ha="left", fontsize=9.5, fontweight="semibold")
    save(fig, name, "Shaded: within 5 %. Water: H2O at 1 g/cm³ by Bragg's rule. Bethe-Bloch drawn from 1 MeV/u.")


@figure("nist_pstar")
def fig_pstar(inp):
    nist_figure(inp, "pstar", "H", 1.00728, "nist_pstar", "Protons against NIST PSTAR (ICRU 49)")


@figure("nist_astar")
def fig_astar(inp):
    nist_figure(inp, "astar", "He", 4.0015, "nist_astar", "Alpha particles against NIST ASTAR (ICRU 49)")


@figure("nist_csda")
def fig_csda(inp):
    import stopping_reference as sr
    fig, axes = plt.subplots(2, 3, figsize=(9.5, 5.4), sharex=True, sharey=True)
    for ax, (mat, target, density) in zip(axes.flat, NIST_MATS):
        ratio_axes(ax, 0.75, 1.25, band=0.05)
        for prog, ion, mass, color in (("pstar", "H", 1.00728, BLUE), ("astar", "He", 4.0015, ORANGE)):
            ref = nist_rows(inp.need(inp.a.nist / f"{prog}_{mat}.txt"))
            lo, hi = ref[0][0] * 1000 / 1.2, ref[-1][0] * 1000 * 1.2
            t = inp.table(ion, target, lo, hi, "srim", density=density)
            if density:
                rho = float(density)
            else:
                rho = float(subprocess_density(inp, target))
            pts = [(r[0] / mass, sr.interpolate(t, r[0] * 1000.0, 3) / (r[4] / rho * 1e4)) for r in ref]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=1.6, label=f"{ion} ({prog.upper()})")
        ax.set_xscale("log")
        ax.set_title(mat)
    for ax in axes[1]:
        ax.set_xlabel("energy per mass (MeV/u)")
    for ax in axes[:, 0]:
        ax.set_ylabel("CSDA range, AMCSET / NIST")
    axes[0, 0].legend(loc="lower right", fontsize=7)
    save(fig, "nist_csda", "Continuous-slowing-down path length from the total (electronic + nuclear) stopping of the bundled tables. Shaded: within 5 %.")


def subprocess_density(inp, target):
    import subprocess
    head = subprocess.run([str(inp.a.exe), "stopping", "--ion", "H", "--target", target, "--from", "1 MeV",
                           "--to", "2 MeV", "--per-decade", "1"], capture_output=True, text=True).stdout
    return head.split("density ")[1].split(" g/cm3")[0]


# =========================================================================================
# 3. Ni ions
# =========================================================================================

@figure("ni_peaks_energy")
def fig_ni_peaks(inp):
    ni = inp.lit()["ni_peaks"]
    fig, axes = plt.subplots(1, 3, figsize=(10.0, 3.7), gridspec_kw={"width_ratios": [1, 1.3, 0.7]})
    for ax, (title, sel) in zip(axes, (("62Ni into Ni", [r for r in ni if r["target"] == "Ni"]),
                                        ("Ni into Fe", [r for r in ni if r["target"] == "Fe"]),
                                        ("Ni into steels", [r for r in ni if r["target"] in ("316", "304")]))):
        sel = sorted(sel, key=lambda r: r["energy_MeV"])
        e = [r["energy_MeV"] for r in sel]
        ax.plot(e, [r["srim_um"] for r in sel], "o-", color=ORANGE, label="AMCSET, SRIM-2013 tables")
        ax.plot(e, [r["fix_um"] for r in sel], "D-", color=BLUE, ms=5, label="AMCSET with nickel_ions")
        sp = [(r["energy_MeV"], r["fix_sputtered_um"]) for r in sel if r["fix_sputtered_um"]]
        ax.plot([x for x, _ in sp], [y for _, y in sp], "D", color=BLUE, mfc=SURF, mec=BLUE, mew=1.2, ms=5,
                label="  … seen after sputtering")
        ax.errorbar(e, [r["measured_um"] for r in sel], yerr=[r["measured_err_um"] or 0 for r in sel],
                    ecolor=INK, elinewidth=1, **meas_style(ms=5.5), label="measured peak")
        ax.set_xlabel("ion energy (MeV)")
        ax.set_title(title)
        ax.set_ylim(0, None)
        if ax is axes[2]:
            ax.set_xticks(e)
            ax.set_xticklabels([f"{r['target']} SS\n{r['energy_MeV']:g} MeV" for r in sel])
            ax.set_xlim(7, 10)
            ax.set_xlabel("")
    axes[0].set_ylabel("peak depth of the implanted Ni (µm)")
    axes[1].legend(loc="upper left", fontsize=7.5)
    save(fig, "ni_peaks_energy", "Measured: Pena et al. 2025 (Ni, 2.0 and 3.5 MeV), Chauhan et al. 2026 (Ni 5 MeV; Fe 2-5 MeV), "
         "Song et al. 2026 Table 2 (Fe 8 MeV after Lin et al., 11.5 MeV after Gurbich et al.; 316 and 304 SS).")


@figure("ni_peaks_ratio")
def fig_ni_ratio(inp):
    ni = inp.lit()["ni_peaks"]
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    vratio_axes(ax)
    for i, r in enumerate(ni):
        y = len(ni) - 1 - i
        m = r["measured_um"]
        ax.plot(r["srim_um"] / m, y, "o", color=ORANGE, ms=6)
        ax.plot(r["fix_um"] / m, y, "D", color=BLUE, ms=5.5)
        if r["fix_sputtered_um"]:
            ax.plot(r["fix_sputtered_um"] / m, y, "D", color=BLUE, mfc=SURF, mec=BLUE, mew=1.2, ms=5)
        if r["measured_err_um"]:
            e = r["measured_err_um"] / m
            ax.plot([1 - e, 1 + e], [y, y], color=INK, lw=1.0)
    ax.set_yticks(range(len(ni)))
    ax.set_yticklabels([f"{r['energy_MeV']:g} MeV → {r['target'] + (' SS' if r['target'] in ('316', '304') else '')}"
                        for r in ni][::-1])
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("AMCSET peak / measured peak")
    ax.legend([Line2D([], [], marker="o", color=ORANGE, ls="none"), Line2D([], [], marker="D", color=BLUE, ls="none"),
               Line2D([], [], marker="D", color=BLUE, mfc=SURF, mec=BLUE, mew=1.2, ls="none"),
               Line2D([], [], color=INK, lw=1)],
              ["SRIM-2013 tables", "nickel_ions", "nickel_ions, after sputtering", "measurement's uncertainty"],
              loc="upper center", bbox_to_anchor=(0.45, -0.13), ncol=2, fontsize=7.5)
    ax.set_title("Ni peak depths: AMCSET over measured")
    save(fig, "ni_peaks_ratio", "Shaded: within 10 %.")


@figure("ni_profiles")
def fig_ni_profiles(inp):
    ni = {r["key"]: r for r in inp.lit()["ni_peaks"]}
    keys = [("Ni_3.5", "3.5 MeV 62Ni → Ni"), ("Ni_5.0", "5 MeV 62Ni → Ni"), ("Fe_5.0", "5 MeV 62Ni → Fe"),
            ("316_8.0", "8 MeV Ni → 316 SS")]
    fig, axes = plt.subplots(1, 4, figsize=(10.5, 3.0), sharey=True)
    for ax, (k, title) in zip(axes, keys):
        r = ni[k]
        x = r["profile_depth_um"]
        for prof, color, label in ((r["profile_srim"], ORANGE, "SRIM-2013 tables"), (r["profile_fix"], BLUE, "nickel_ions")):
            mx = max(prof)
            ax.plot(x, [v / mx for v in prof], color=color, lw=1.6, label=label)
        m, e = r["measured_um"], r["measured_err_um"] or 0
        if e:
            ax.axvspan(m - e, m + e, color=BAND, lw=0, zorder=0)
        ax.axvline(m, color=INK, lw=1.2, label="measured peak")
        ax.set_title(title)
        ax.set_xlabel("depth (µm)")
        ax.set_xlim(0, max(x))
    axes[0].set_ylabel("implanted Ni (per maximum)")
    axes[0].legend(loc="upper right", fontsize=7)
    save(fig, "ni_profiles", "Stopped-ion profiles of AMCSET (quick calculation, 4 × 20 000 ions) and the measured peak (shaded: its uncertainty).")


@figure("ni_sputtering")
def fig_ni_sputter(inp):
    ni = [r for r in inp.lit()["ni_peaks"] if r["sputtered_nm"]]
    fig, ax = plt.subplots(figsize=(6.0, 3.2))
    labels = [f"{r['energy_MeV']:g} MeV → {r['target']}" for r in ni]
    ax.barh(range(len(ni)), [r["sputtered_nm"] for r in ni], color=BLUE, height=0.6)
    for i, r in enumerate(ni):
        ax.text(r["sputtered_nm"] + 2, i, f"Y = {r['sputtering_yield']:.2f}, {r['fluence_cm2']:.2g} /cm²",
                va="center", fontsize=7.5, color=INK2)
    ax.set_yticks(range(len(ni)))
    ax.set_yticklabels(labels)
    ax.invert_yaxis()
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 160)
    ax.set_xlabel("depth sputtered away during the implantation, d = Y × fluence / N (nm)")
    ax.set_title("How much surface the Ni implantations removed")
    save(fig, "ni_sputtering", "Y: AMCSET full-cascade sputtering yield with nickel_ions. The measured peaks are not corrected for it.")


# =========================================================================================
# 4. Ranges
# =========================================================================================

GROUP_LABEL = {"behar_si": "a-Si (Behar, Tables 3-4)", "behar_sio2": "SiO2 (Behar, Table 5)",
               "behar_light": "H, Li, B in a-Si (Behar, 6-9)", "oosterhoff": "B, P in Si (Oosterhoff)",
               "behar_c": "C films (Behar, 11-13)", "behar_b": "B films (Behar, 15)",
               "behar_sic": "SiC films (Behar, 16)", "zhang_sic": "SiC, peaks (Zhang 2009)"}
HEAVY = 55


@figure("ranges_parity")
def fig_ranges_parity(inp):
    rg = inp.ranges()
    sets = [("Si and SiO2", ("behar_si", "behar_sio2", "behar_light", "oosterhoff"), BLUE, "o"),
            ("C and B films", ("behar_c", "behar_b"), ORANGE, "s"),
            ("SiC", ("behar_sic", "zhang_sic"), AQUA, "^")]
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    xs = [5, 6000]
    ax.fill_between(xs, [x * 0.9 for x in xs], [x * 1.1 for x in xs], color=BAND, lw=0, zorder=0)
    ax.plot(xs, xs, color=AXIS, lw=0.9)
    for label, groups, color, mk in sets:
        sel = [r for r in rg if r["group"] in groups]
        ax.plot([float(r["measured_nm"]) for r in sel], [float(r["amcset_nm"]) for r in sel], mk, color=color,
                ms=4.5, label=f"{label} ({len(sel)})")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xs)
    ax.set_ylim(*xs)
    ax.set_xlabel("measured range (nm)")
    ax.set_ylabel("AMCSET range (nm)")
    ax.legend(loc="upper left")
    ax.set_title("350 measured ranges against AMCSET")
    save(fig, "ranges_parity", "Mean projected ranges (Zhang et al.: peak depths). Shaded: within 10 %.")


@figure("ranges_energy")
def fig_ranges_energy(inp):
    rg = inp.ranges()
    import stopping_reference as sr
    fig, axes = plt.subplots(2, 4, figsize=(11.0, 5.2), sharey=True)
    for ax, g in zip(axes.flat, GROUP_LABEL):
        sel = [r for r in rg if r["group"] == g]
        ratio_axes(ax, 0.5, 1.6)
        for heavy, color, mk in ((False, BLUE, "o"), (True, ORANGE, "s")):
            s = [r for r in sel if (sr.Z[r["ion"]] >= HEAVY) == heavy]
            ax.plot([float(r["energy_keV"]) for r in s], [float(r["ratio"]) for r in s], mk, color=color, ms=4)
        ax.set_xscale("log")
        ax.set_title(GROUP_LABEL[g], fontsize=8)
    for ax in axes[1]:
        ax.set_xlabel("energy (keV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("AMCSET / measured")
    axes[0, 0].legend([Line2D([], [], marker="o", color=BLUE, ls="none"), Line2D([], [], marker="s", color=ORANGE, ls="none")],
                      [f"Z1 < {HEAVY}", f"Z1 ≥ {HEAVY}"], loc="lower left", fontsize=7)
    save(fig, "ranges_energy", "Ratio of AMCSET's mean projected range to the measured one, by data set. Shaded: within 10 %.")


@figure("ranges_z1")
def fig_ranges_z1(inp):
    rg = inp.ranges()
    import stopping_reference as sr
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    ratio_axes(ax, 0.65, 1.3)
    for g, color, label in (("behar_si", BLUE, "a-Si"), ("behar_sio2", AQUA, "SiO2"), ("behar_c", ORANGE, "C"),
                            ("behar_b", MAGENTA, "B")):
        by = {}
        for r in rg:
            if r["group"] == g:
                by.setdefault(r["ion"], []).append(float(r["ratio"]))
        pts = sorted((sr.Z[i], statistics.median(v)) for i, v in by.items())
        ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=color, ms=4.5, lw=1.3, label=label)
        ax.text(pts[-1][0] + 1.5, pts[-1][1], label, color=INK2, fontsize=7.5, va="center")
    ax.set_xlabel("ion Z1")
    ax.set_ylabel("median AMCSET / measured R_p")
    ax.set_xlim(15, 92)
    ax.legend(loc="lower left", ncol=4, fontsize=7.5)
    ax.set_title("Range ratio by ion: the light-target deficit grows with Z1")
    save(fig, "ranges_z1", "Behar et al. 1995; median over each ion's energies (10-400 keV). Shaded: within 10 %.")


@figure("ranges_straggling")
def fig_ranges_straggling(inp):
    rg = inp.ranges()
    groups = [g for g in GROUP_LABEL if g != "zhang_sic"]
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    vratio_axes(ax)
    for i, g in enumerate(groups):
        v = [float(r["straggling_ratio"]) for r in rg if r["group"] == g and r["straggling_ratio"]]
        y = len(groups) - 1 - i
        ax.plot(v, [y] * len(v), "o", color=BLUE, alpha=0.45, mew=0, ms=4)
        ax.plot(statistics.median(v), y, "|", color=INK, ms=14, mew=2)
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([GROUP_LABEL[g] for g in groups][::-1], fontsize=7.5)
    ax.grid(axis="y", visible=False)
    ax.set_xscale("log")
    ax.set_xticks([0.3, 0.5, 0.7, 1, 1.4, 2])
    ax.set_xticklabels(["0.3", "0.5", "0.7", "1", "1.4", "2"])
    ax.set_xlabel("AMCSET / measured range straggling (bar: median)")
    ax.set_title("Range straggling")
    save(fig, "ranges_straggling", "Shaded: within 10 %. Measured straggling includes the depth resolution of RBS and the profiles' deconvolution.")


@figure("oosterhoff")
def fig_oosterhoff(inp):
    rg = [r for r in inp.ranges() if r["group"] == "oosterhoff"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.0, 3.4))
    for ion, color in (("B", BLUE), ("P", ORANGE)):
        s = sorted((r for r in rg if r["ion"] == ion), key=lambda r: float(r["energy_keV"]))
        e = [float(r["energy_keV"]) for r in s]
        a.plot(e, [float(r["amcset_nm"]) for r in s], "-", color=color, label=f"{ion}, AMCSET")
        a.plot(e, [float(r["measured_nm"]) for r in s], "o", color=color, label=f"{ion}, SIMS")
        b.plot(e, [float(r["amcset_straggling_nm"]) for r in s], "-", color=color)
        b.plot(e, [float(r["measured_straggling_nm"]) for r in s], "o", color=color)
    a.set_xlabel("energy (keV)")
    a.set_ylabel("projected range R_p (nm)")
    a.legend(loc="upper left", fontsize=7.5)
    a.set_title("Projected range")
    b.set_xlabel("energy (keV)")
    b.set_ylabel("straggling ΔR_p (nm)")
    b.set_title("Straggling")
    save(fig, "oosterhoff", "B and P in (100) Si at 7° tilt, as-implanted SIMS (Oosterhoff 1988, Table 1); AMCSET with an amorphous target.")


@figure("light_ions_si")
def fig_light(inp):
    rg = [r for r in inp.ranges() if r["group"] == "behar_light"]
    fig, ax = plt.subplots(figsize=(6.0, 3.8))
    for ion, color in (("H", BLUE), ("Li", ORANGE), ("B", AQUA)):
        s = sorted((r for r in rg if r["ion"] == ion), key=lambda r: float(r["energy_keV"]))
        e = [float(r["energy_keV"]) for r in s]
        ax.plot(e, [float(r["amcset_nm"]) for r in s], "-", color=color)
        ax.plot(e, [float(r["measured_nm"]) for r in s], "o", color=color, label=ion)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("energy (keV)")
    ax.set_ylabel("projected range R_p (nm)")
    ax.legend(loc="upper left", title="points: measured, lines: AMCSET", title_fontsize=7)
    ax.set_title("H, Li and B in amorphous Si (Behar et al. 1995, Tables 6-9)")
    save(fig, "light_ions_si", "Nuclear-reaction depth profiling (15N resonance for H; neutron depth profiling and (p,γ) for Li and B).")


@figure("zhang_sic")
def fig_zhang(inp):
    import range_data as rd
    rg = [r for r in inp.ranges() if r["group"] == "zhang_sic"]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.5))
    for ax, ion in zip(axes, ("Pb", "Au")):
        s = sorted((r for r in rg if r["ion"] == ion), key=lambda r: float(r["energy_keV"]))
        e = [float(r["energy_keV"]) for r in s]
        ax.plot(e, [rd.ZHANG_SRIM_PEAK[(ion, int(float(r["energy_keV"])))] for r in s], "s", color=ORANGE,
                mfc=SURF, mec=ORANGE, mew=1.3, label="their SRIM-2008")
        ax.plot(e, [float(r["amcset_nm"]) for r in s], "-", color=BLUE, label="AMCSET")
        ax.plot(e, [float(r["measured_nm"]) for r in s], **meas_style(ms=5.5), label="measured (RBS, SIMS)")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("energy (keV)")
        ax.set_title(f"{ion} in SiC: peak depth")
    axes[0].set_ylabel("peak depth (nm)")
    axes[1].legend(loc="upper left", fontsize=7.5)
    save(fig, "zhang_sic", "Zhang et al., J. Appl. Phys. 105, 104901 (2009), Tables I-III; 2 MeV Au at 60°, the rest near normal incidence.")


# =========================================================================================
# 5. Sputtering
# =========================================================================================

def sputter_panels(inp, pairs, name, title):
    sp = inp.sputter()
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 5.6))
    for ax, (ion, target) in zip(axes.flat, pairs):
        s = sorted((r for r in sp if r["ion"] == ion and r["target"] == target), key=lambda r: float(r["energy_keV"]))
        e = [float(r["energy_keV"]) for r in s]
        a = [(float(r["energy_keV"]), float(r["amcset"])) for r in s if float(r["amcset"]) > 0]
        ax.plot([x for x, _ in a], [y for _, y in a], "-", color=BLUE, label="AMCSET (full cascades)")
        ax.plot([r for r in e], [float(r["measured"]) for r in s], **meas_style(), label="measured")
        out = [r for r in s if r["flag"]]
        ax.plot([float(r["energy_keV"]) for r in out], [float(r["measured"]) for r in out], "x", color=RED, mec=RED,
                ms=7, mew=1.4, label="low outlier (left out of the ratios)")
        zero = [float(r["energy_keV"]) for r in s if float(r["amcset"]) == 0]
        if zero:  # on the bottom of the axes: no yield to draw on a log axis
            ax.plot(zero, [0.03] * len(zero), "v", color=BLUE, ms=5, transform=ax.get_xaxis_transform())
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{ion} on {target}")
    for ax in axes[1]:
        ax.set_xlabel("ion energy (keV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("sputtering yield (atoms/ion)")
    axes[0, 0].legend(loc="lower right", fontsize=7)
    fig.suptitle(title, x=0.0, y=1.03, ha="left", fontsize=9.5, fontweight="semibold")
    save(fig, name, "Measured: points read off the figures of Yamamura et al. (NIFS-DATA-23, -31), 20-30 % reading error. "
         "Triangles on the axis: energies where AMCSET sputters nothing.")


@figure("sputter_argon")
def fig_sputter_ar(inp):
    sputter_panels(inp, [("Ar", t) for t in ("Si", "Fe", "Cu", "Ag", "W", "Au")], "sputter_argon",
                   "Sputtering yields of Ar ions at normal incidence")


@figure("sputter_light")
def fig_sputter_light(inp):
    sputter_panels(inp, [(i, t) for t in ("Fe", "W") for i in ("H", "D", "He")], "sputter_light",
                   "Sputtering yields of H, D and He ions at normal incidence")


@figure("sputter_summary")
def fig_sputter_summary(inp):
    sp = [r for r in inp.sputter() if not r["flag"]]
    pairs = []
    for r in sp:
        if (r["ion"], r["target"]) not in pairs:
            pairs.append((r["ion"], r["target"]))
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    vratio_axes(ax)
    for i, (ion, t) in enumerate(pairs):
        y = len(pairs) - 1 - i
        hi = [float(r["amcset"]) / float(r["measured"]) for r in sp if r["ion"] == ion and r["target"] == t
              and float(r["energy_keV"]) >= 0.5 and float(r["amcset"]) > 0]
        lo = [float(r["amcset"]) / float(r["measured"]) for r in sp if r["ion"] == ion and r["target"] == t
              and float(r["energy_keV"]) < 0.5 and float(r["amcset"]) > 0]
        ax.plot(lo, [y] * len(lo), "o", color=ORANGE, ms=4, alpha=0.7, mew=0)
        ax.plot(hi, [y] * len(hi), "o", color=BLUE, ms=4.5, alpha=0.8, mew=0)
        if hi:
            ax.plot(statistics.median(hi), y, "|", color=INK, ms=14, mew=2)
    ax.set_yticks(range(len(pairs)))
    ax.set_yticklabels([f"{i} on {t}" for i, t in pairs][::-1])
    ax.grid(axis="y", visible=False)
    ax.set_xscale("log")
    ax.set_xticks([0.1, 0.2, 0.5, 1, 2, 5])
    ax.set_xticklabels(["0.1", "0.2", "0.5", "1", "2", "5"])
    ax.set_xlabel("AMCSET / measured sputtering yield")
    ax.legend([Line2D([], [], marker="o", color=BLUE, ls="none"), Line2D([], [], marker="o", color=ORANGE, ls="none"),
               Line2D([], [], marker="|", color=INK, ls="none", ms=10, mew=2)],
              ["E ≥ 0.5 keV", "E < 0.5 keV", "median, E ≥ 0.5 keV"], loc="lower left", fontsize=7.5)
    ax.set_title("Sputtering: each ion-target pair has its own factor")
    save(fig, "sputter_summary", "Shaded: within 10 %. Zero yields (near the light-ion thresholds) cannot be drawn on a log axis and are left out here.")


# =========================================================================================
# 6. Damage
# =========================================================================================

@figure("damage_functions")
def fig_damage(inp):
    rows = read_csv(inp.need(inp.a.damage / "eb0" / "damage_functions.csv"))
    els = []
    for r in rows:
        if r["element"] not in els:
            els.append(r["element"])
    fig, axes = plt.subplots(2, 3, figsize=(10.0, 5.4), sharex=True, sharey=True)
    for ax, el in zip(axes.flat, els):
        s = sorted((r for r in rows if r["element"] == el), key=lambda r: float(r["E_keV"]))
        e = [float(r["E_keV"]) for r in s]
        nrt = [float(r["nrt"]) for r in s]
        for key, color, label in (("trim_vac", ORANGE, "TRIM full-cascade vacancies"), ("kp", AQUA, "Kinchin-Pease"),
                                  ("arc", BLUE, "arc-dpa")):
            ax.plot(e, [float(r[key]) / n for r, n in zip(s, nrt)], "o-", color=color, ms=3.8, lw=1.4, label=label)
        ax.axhline(1.0, color=AXIS, lw=0.9)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(f"{el} (Ed {float(s[0]['Ed_eV']):g} eV)")
    for ax in axes[1]:
        ax.set_xlabel("PKA energy (keV)")
    for ax in axes[:, 0]:
        ax.set_ylabel("per PKA, / NRT")
    axes[0, 0].legend(loc="lower left", fontsize=7)
    ax = axes[0, 0]
    ax.set_yticks([0.1, 0.2, 0.5, 1, 2, 3])
    ax.set_yticklabels(["0.1", "0.2", "0.5", "1", "2", "3"])
    save(fig, "damage_functions", "PKAs of one energy deep in pure targets, full cascades, Eb = 0 (pka_damage_functions.py). "
         "arc-dpa is a fit to MD cascades (Nordlund et al. 2018); NRT = 1.")


@figure("mitsi_parity")
def fig_mitsi(inp):
    txt = inp.need(inp.a.damage / "mitsi2023" / "comparison.txt").read_text(encoding="utf-8")
    pts = [(float(m.group(3)), float(m.group(4)), m.group(1), m.group(2)) for m in
           re.finditer(r"^(\w+)\s+(\w+ \d+ MeV)\s+([\d.]+) \+-\s*[\d.]+\s+([\d.]+)", txt, re.M)]
    fig, ax = plt.subplots(figsize=(4.8, 4.4))
    ax.plot([0.1, 0.8], [0.1, 0.8], color=AXIS, lw=0.9)
    beams = []
    for _, _, _, b in pts:
        if b not in beams:
            beams.append(b)
    for k, b in enumerate(beams):
        s = [p for p in pts if p[3] == b]
        ax.plot([p[1] for p in s], [p[0] for p in s], "o", color=C[k], ms=5, label=b)
    ax.set_xlabel("arc-dpa / NRT, Mitsi et al. 2023 (SRIM-2013 post-processing)")
    ax.set_ylabel("arc-dpa / NRT, AMCSET")
    ax.legend(loc="upper left", fontsize=7.5)
    ax.set_title("arc-dpa of whole ion tracks: 30 cases")
    save(fig, "mitsi_parity", "Targets Fe, Ni, Cu, Pd, W, Pt; rms difference 0.43 %.")


# =========================================================================================
# 7. Crystals
# =========================================================================================

def crystal_text(inp):
    return inp.need(inp.a.crystal / "results.txt").read_text(encoding="utf-8")


@figure("crystal_dip")
def fig_dip(inp):
    txt = crystal_text(inp)
    pts = [(float(a), float(b)) for a, b in re.findall(r"tilt ([\d.]+) deg\s+[\d.]+\s+chi =\s+([\d.]+)", txt)]
    fig, ax = plt.subplots(figsize=(5.6, 3.4))
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=BLUE, label="AMCSET, crystal transport")
    ax.axhline(1.0, color=AXIS, lw=0.9)
    ax.plot(0.0, 0.03, **meas_style(marker="D", ms=6), label="measured minimum yield (0.03)")
    ax.errorbar(0.95, 0.51, xerr=0.0, **meas_style(marker="s", ms=6), label="measured half-angle (0.95°)")
    ax.set_xlabel("tilt from ⟨110⟩ (degrees)")
    ax.set_ylabel("normalized yield χ")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.set_title("Channeling dip: 1 MeV He along ⟨110⟩ of Ge, 293 K")
    save(fig, "crystal_dip", "Atoms displaced by the ion at 20-120 nm relative to amorphous Ge (calibration.py); measurement: Picraux et al., as quoted by Gemmell 1974.")


@figure("crystal_boron")
def fig_boron(inp):
    txt = crystal_text(inp)
    levels = [1e5, 1e4, 1e3, 1e2]
    rows = {}
    for line in txt.split("0.5 keV B into")[1].splitlines():
        m = re.match(r"\s+(.+?)\s{2,}([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", line)
        if m:
            rows[m.group(1).strip()] = [float(m.group(i)) for i in range(2, 6)]
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    series = [("measured, 1e14 cm-2", INK, "o", "measured (10¹⁴ cm⁻²)"), ("trim2013 (amorphous)", ORANGE, "s", "amorphous (trim2013)"),
              ("crystal 0 deg, fraction 0", RAMP[250], "^", "crystal, division 0"),
              ("crystal 0 deg, fraction 0.5", RAMP[450], "^", "crystal, division 0.5 (default)"),
              ("crystal 0 deg, fraction 1", RAMP[650], "^", "crystal, division 1")]
    for key, color, mk, label in series:
        if key in rows:
            ax.plot(rows[key], levels, mk + "-", color=color, ms=5, lw=1.4, label=label)
    ax.set_yscale("log")
    ax.set_xlabel("depth (nm)")
    ax.set_ylabel("concentration / dose (1/cm)")
    ax.legend(loc="upper right", fontsize=7.5)
    ax.set_title("0.5 keV B into (100) Si: depth at each concentration level")
    save(fig, "crystal_boron", "Schroer et al. 2000, Fig. 1 (read off the figure). 'Division': the share of electronic stopping lost in collisions.")


@figure("as_channeling")
def fig_as(inp):
    rows = {r["label"]: r for r in inp.lit()["as_si_channeling"]}
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6), sharey=True)
    for ax, tilt in zip(axes, (0, 10)):
        for label, color in (("amorphous Si", ORANGE),
                             ("crystal, tilt 0" if tilt == 0 else "crystal, tilt 10 deg", BLUE),
                             (f"1 nm SiO2 + crystal, {tilt}", AQUA)):
            r = rows[label]
            pts = [(x, v) for x, v in zip(r["profile_depth_nm"], r["profile_over_peak"]) if v > 0]
            ax.plot([p[0] for p in pts], [p[1] for p in pts], color=color, lw=1.4,
                    label=label.replace(", tilt 0", "").replace(", tilt 10 deg", "").replace(f", {tilt}", ""))
        m = rows["crystal, tilt 0" if tilt == 0 else "crystal, tilt 10 deg"]
        ax.plot([12, 100, 200], [1.0, m["measured_c100"], m["measured_c200"]], **meas_style(ms=6), label="measured (Tasch et al.)")
        ax.set_yscale("log")
        ax.set_ylim(1e-6, 2)
        ax.set_xlim(0, 300)
        ax.set_xlabel("depth (nm)")
        ax.set_title(f"tilt {tilt}°")
    axes[0].set_ylabel("As concentration / peak")
    axes[0].legend(loc="upper right", fontsize=7.5)
    fig.suptitle("15 keV As into (100) Si, 10¹⁴ cm⁻²: channeling tails", x=0.0, y=1.03, ha="left", fontsize=9.5, fontweight="semibold")
    save(fig, "as_channeling", "Measured: read off Tasch et al. 1994, Fig. 3 (about a factor 2). AMCSET: 200 000 ions, perfect crystal, default stopping division; "
         "no ion of the amorphous calculation stops beyond about 50 nm.")


# =========================================================================================
# 8. Ionization
# =========================================================================================

@figure("ionization_yield")
def fig_quenching(inp):
    q = inp.lit()["quenching"]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.6))
    for ax, target in zip(axes, ("Ge", "Si")):
        s = sorted((r for r in q if r["target"] == target), key=lambda r: r["energy_keV"])
        e = [r["energy_keV"] for r in s]
        ax.plot(e, [r["amcset"] for r in s], "o-", color=BLUE, ms=4, label="AMCSET, SRIM tables")
        ax.plot(e, [r["amcset_lss"] for r in s], "s-", color=ORANGE, ms=4, label="AMCSET, LSS stopping")
        ax.plot(e, [r["lindhard"] for r in s], "-", color=AQUA, lw=1.4, label="Lindhard")
        for src, mk in (("Bonhomme 2022", "o"), ("Chavarria 2016", "o"), ("Izraelevitch 2017", "s")):
            m = [r for r in s if r["source"] == src]
            if m:
                ax.errorbar([r["energy_keV"] for r in m], [r["measured"] for r in m], yerr=[r["measured_err"] for r in m],
                            ecolor=INK, elinewidth=1, **meas_style(marker=mk, ms=5), label=f"measured ({src})")
        ax.set_xscale("log")
        ax.set_xlabel("recoil energy (keV)")
        ax.set_title(f"{target} recoils in {target}")
        ax.set_ylim(0, 0.34 if target == "Ge" else None)
    axes[0].set_ylabel("ionization / recoil energy")
    axes[0].legend(loc="lower right", fontsize=7)
    axes[1].legend(loc="upper left", fontsize=7)
    save(fig, "ionization_yield", "PKAs started deep in a pure target, full cascades (500 PKAs × 4 seeds per point).")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    b = ROOT / "build"
    ap.add_argument("--exe", type=Path, default=b / "Release" / ("amcset.exe" if sys.platform == "win32" else "amcset"))
    ap.add_argument("--literature", type=Path, default=b / "validation" / "literature")
    ap.add_argument("--stopping", type=Path, default=b / "validation" / "stopping")
    ap.add_argument("--scan", type=Path, default=b / "validation" / "table_scan")
    ap.add_argument("--nist", type=Path, default=b / "validation" / "sources" / "data")
    ap.add_argument("--trim", type=Path, default=b / "trim_validation")
    ap.add_argument("--damage", type=Path, default=b / "damage_models")
    ap.add_argument("--crystal", type=Path, default=b / "crystal_calibration")
    ap.add_argument("--only", default="")
    ap.add_argument("--png", type=Path, help="also write PNG previews into this folder")
    a = ap.parse_args()
    global PNG
    PNG = a.png
    if not a.exe.exists():
        a.exe = None
    inp = Inputs(a)
    names = a.only.split(",") if a.only else list(FIGURES)
    missing = []
    for n in names:
        try:
            FIGURES[n](inp)
        except FileNotFoundError as e:
            missing.append(f"{n}: missing {e}")
            plt.close("all")
    for m in missing:
        print("skipped", m)


if __name__ == "__main__":
    main()

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
"""Figures of docs/GPU.md from the outputs of compare_backends.py and speed.py.

Reads (default --data build/gpu_compare):
  double/analysis.json     OpenCL double against the CPU (compare_backends.py)
  mixed/analysis.json      mixed precision before the precision fix
  mixed3/analysis.json     mixed precision after it
  more/<case>/<backend>/   larger checks (10 runs each, summary.json)
  refcheck/big_*/          Fe 50 keV quick damage, 10 x 1,000,000 ions
  speed/speed.json         speed.py
  live/<N>/                speed.py --live hip:mixed (the first N ions of one run)
and writes SVG figures and two GIF animations to docs/gpu/figures/. The
live-plot animation is data (each frame is what the GUI's plots show after
N ions); the batches-against-ring animation is an illustration: its block
lengths are drawn from a distribution with the measured spread.
Needs matplotlib and numpy (only this script; the evidence comes from the
standard-library tools).

  python tools/gpu_compare/make_figures.py [--data build/gpu_compare] [--out docs/gpu/figures]
"""

import argparse
import glob
import json
import math
import os
import statistics
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.transforms  # noqa: E402,F401
import numpy as np  # noqa: E402
from matplotlib import animation  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_backends import quantities  # noqa: E402

INK = "#0f172a"
MUTED = "#64748b"
GRID = "#e2e8f0"
COLORS = {"cpu": "#2563eb", "opencl:double": "#f59e0b", "opencl:mixed": "#ea580c",
          "hip:double": "#a78bfa", "hip:mixed": "#7c3aed", "reference:mixed": "#0d9488"}
LABELS = {"cpu": "CPU (16 threads)", "opencl:double": "OpenCL double", "opencl:mixed": "OpenCL mixed",
          "hip:double": "HIP double", "hip:mixed": "HIP mixed"}
CASE_LABELS = {"fe50keV_full": "Fe 50 keV, full", "fe50keV_kp": "Fe 50 keV, quick",
               "fe1MeV_full": "Fe 1 MeV, full", "ar30keV_SiO2_Si_full": "Ar 30 keV SiO2/Si",
               "h100keV_Si_full": "H 100 keV Si", "au200keV_Si_full": "Au 200 keV Si",
               "ar2keV_Cu_monolayer": "Ar 2 keV Cu, monolayer"}


def style():
    plt.rcParams.update({
        "font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 10,
        "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": GRID,
        "grid.linewidth": 0.8, "axes.axisbelow": True, "legend.frameon": False, "svg.fonttype": "none",
        "figure.facecolor": "white", "axes.facecolor": "white", "lines.linewidth": 2,
    })


def t_pdf(x, dof):
    return math.exp(math.lgamma((dof + 1) / 2) - math.lgamma(dof / 2)
                    - (dof + 1) / 2 * math.log1p(x * x / dof)) / math.sqrt(dof * math.pi)


def zs(path, backend):
    d = json.load(open(path))
    out = []
    for case, v in d.items():
        if case.startswith("_"):
            continue
        for k, q in v["quantities"].items():
            # A quantity constant in every run (both errors 0) carries no information.
            if backend in q and math.isfinite(q[backend][2]) and math.hypot(q[backend][1], q["cpu"][1]) > 0:
                out.append((q[backend][2], case, k))
    return out


def save(fig, out, name):
    p = os.path.join(out, name)
    fig.savefig(p, bbox_inches="tight", metadata={"Date": None})  # no date: redrawn files differ only by data
    plt.close(fig)
    with open(p, "rb") as f:
        data = f.read()
    with open(p, "wb") as f:  # LF line endings (the repository's), also on Windows
        f.write(data.replace(b"\r\n", b"\n"))
    print("wrote", p)


def fig_z_distribution(data, out):
    panels = [("double", "opencl:double", "OpenCL double"),
              ("mixed", "opencl:mixed", "OpenCL mixed, before the fix"),
              ("mixed3", "opencl:mixed", "OpenCL mixed, after"),
              ("mixed", "hip:mixed", "HIP mixed, before the fix"),
              ("mixed3", "hip:mixed", "HIP mixed, after")]
    fig, axes = plt.subplots(1, len(panels), figsize=(14, 3.2), sharey=True)
    x = np.linspace(-4.5, 4.5, 300)
    for ax, (run, b, title) in zip(axes, panels):
        p = os.path.join(data, run, "analysis.json")
        if not os.path.exists(p):
            ax.set_visible(False)
            continue
        z = np.array([t[0] for t in zs(p, b)])
        bins = np.arange(-4.5, 4.75, 0.5)
        ax.hist(np.clip(z, -4.49, 4.49), bins=bins, density=True, color=COLORS[b], alpha=0.75,
                edgecolor="white", linewidth=1)
        ax.plot(x, [t_pdf(v, 10) for v in x], color=INK, linewidth=1.5, label="t, 10 dof (6 vs 6 runs)")
        rms = math.sqrt(float(np.mean(z * z)))
        ax.set_title(f"{title}\n{len(z)} quantities, rms z {rms:.2f}", fontsize=10, color=INK)
        ax.set_xlabel("z = (GPU − CPU) / combined error")
        for v in (-3, 3):
            ax.axvline(v, color=MUTED, linewidth=0.8, linestyle=":")
    axes[0].set_ylabel("density")
    axes[0].legend(loc="upper left", fontsize=8)
    fig.suptitle("GPU against CPU, every quantity of summary.json in seven cases (expected rms 1.12)",
                 color=INK, fontsize=11, y=1.06)
    save(fig, out, "z_distribution.svg")


def fig_z_by_case(data, out):
    p = os.path.join(data, "mixed3", "analysis.json")
    if not os.path.exists(p):
        return
    cases = [c for c in json.load(open(p)) if not c.startswith("_")]
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.axhspan(-2, 2, color="#f1f5f9", zorder=0)
    for i, b in enumerate(["opencl:mixed", "hip:mixed"]):
        rows = zs(p, b)
        for j, c in enumerate(cases):
            z = [t[0] for t in rows if t[1] == c]
            xs = j + (i - 0.5) * 0.3 + (np.random.default_rng(j * 7 + i).random(len(z)) - 0.5) * 0.18
            ax.scatter(xs, z, s=14, color=COLORS[b], alpha=0.8, linewidths=0,
                       label=LABELS[b] if j == 0 else None)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels([CASE_LABELS.get(c, c) for c in cases], rotation=15, ha="right")
    ax.set_ylabel("z per quantity")
    ax.set_ylim(-4.6, 4.6)
    ax.legend(loc="lower right", ncols=2, bbox_to_anchor=(1.0, 1.0))
    ax.set_title("Mixed precision after the fix, by case (6 runs each; grey band |z| < 2)", color=INK, fontsize=11,
                 loc="left")
    save(fig, out, "z_by_case.svg")


def runs_of(pattern):
    return [quantities(json.load(open(f))) for f in sorted(glob.glob(pattern))]


def mean_se(vals):
    return statistics.mean(vals), statistics.stdev(vals) / math.sqrt(len(vals))


def fig_checks(data, out):
    """The deviations that the larger checks decided: the precision fix and two scatter cases."""
    panels = []
    # 1. H 100 keV radial range: before, after (6 runs), and 10 x 500k.
    rows = []
    for label, run, b in [("CPU", "mixed", "cpu"), ("OpenCL mixed\nbefore", "mixed", "opencl_mixed"),
                          ("HIP mixed\nbefore", "mixed", "hip_mixed"), ("OpenCL mixed\nafter", "mixed3", "opencl_mixed"),
                          ("HIP mixed\nafter", "mixed3", "hip_mixed")]:
        v = runs_of(os.path.join(data, run, "h100keV_Si_full", b, "*", "summary.json"))
        if v:
            rows.append((label, *mean_se([r["range.radial_A"] for r in v]), b))
    big = []
    for label, b in [("CPU", "cpu"), ("OpenCL\nmixed", "opencl_mixed"), ("HIP\nmixed", "hip_mixed")]:
        v = runs_of(os.path.join(data, "more", "h100keV_Si_full", b, "*", "summary.json"))
        if v:
            big.append((label, *mean_se([r["range.radial_A"] for r in v]), b))
    panels.append(("H 100 keV in Si: mean radial range (Å)", "6 runs × 50,000 ions", rows,
                   "10 runs × 500,000 ions, after", big))
    # 2. Fe 50 keV quick: ionization by recoils, 6 runs x 100k, then 10 x 1M.
    rows = []
    for label, run, b in [("CPU", "double", "cpu"), ("OpenCL\ndouble", "double", "opencl_double"),
                          ("OpenCL\nmixed", "mixed3", "opencl_mixed"), ("HIP\nmixed", "mixed3", "hip_mixed")]:
        v = runs_of(os.path.join(data, run, "fe50keV_kp", b, "*", "summary.json"))
        if v:
            rows.append((label, *mean_se([r["energy_percent.ionization_recoils"] for r in v]), b))
    big = []
    for label, b in [("CPU", "cpu"), ("OpenCL\ndouble", "opencl_double"), ("OpenCL\nmixed", "opencl_mixed"),
                     ("HIP\nmixed", "hip_mixed")]:
        v = runs_of(os.path.join(data, "refcheck", f"big_{b}_*", "summary.json"))
        if v:
            big.append((label, *mean_se([r["energy_percent.ionization_recoils"] for r in v]), b))
    panels.append(("Fe 50 keV quick damage: energy to ionization by recoils (%)", "6 runs × 100,000 ions", rows,
                   "10 runs × 1,000,000 ions", big))
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2))
    for ax, (title, l1, rows, l2, big) in zip(axes, panels):
        n1 = len(rows)
        for i, (label, m, se, b) in enumerate(rows):
            key = b.replace("_", ":")
            ax.errorbar(i, m, yerr=se, fmt="o", color=COLORS.get(key, INK), capsize=4, markersize=6)
        for i, (label, m, se, b) in enumerate(big):
            key = b.replace("_", ":")
            ax.errorbar(n1 + 0.6 + i, m, yerr=se, fmt="s", color=COLORS.get(key, INK), capsize=4, markersize=6)
        if rows:
            ax.axhline(rows[0][1], color=COLORS["cpu"], linewidth=0.8, linestyle=":")
        if big:
            ax.axhline(big[0][1], color=COLORS["cpu"], linewidth=0.8, linestyle="--")
        ax.set_xticks([*range(n1), *[n1 + 0.6 + i for i in range(len(big))]])
        ax.set_xticklabels([r[0] for r in rows] + [r[0] for r in big], fontsize=8.5)
        ax.axvline(n1 - 0.2, color=GRID, linewidth=1.5)
        ax.set_title(title, color=INK, fontsize=10.5, pad=10)
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax + 0.12 * (ymax - ymin))
        tr = matplotlib.transforms.blended_transform_factory(ax.transData, ax.transAxes)
        ax.text((n1 - 1) / 2, 0.97, l1, ha="center", va="top", color=MUTED, fontsize=9, transform=tr)
        ax.text(n1 + 0.6 + (len(big) - 1) / 2, 0.97, l2, ha="center", va="top", color=MUTED, fontsize=9,
                transform=tr)
    fig.suptitle("Deciding the largest deviations with more ions (error bars: standard error of the mean)",
                 color=INK, fontsize=11, y=1.04)
    fig.tight_layout()
    save(fig, out, "checks.svg")


def fig_speed(data, out):
    p = os.path.join(data, "speed", "speed.json")
    if not os.path.exists(p):
        return
    rows = json.load(open(p))
    groups = []
    for r in rows:
        key = (r["case"], r["ions"])
        if key not in groups:
            groups.append(key)
    backends = ["cpu", "opencl:double", "opencl:mixed", "hip:double", "hip:mixed"]
    fig, axes = plt.subplots(1, len(groups), figsize=(2.3 * len(groups), 4.2))
    for ax, key in zip(np.atleast_1d(axes), groups):
        cpu = [r for r in rows if (r["case"], r["ions"]) == key and r["backend"] == "cpu"]
        base = cpu[0]["wall_seconds"] if cpu else None
        xs, labels = [], []
        for i, b in enumerate(backends):
            r = [r for r in rows if (r["case"], r["ions"]) == key and r["backend"] == b]
            if not r:
                continue
            r = r[0]
            wall = r["wall_seconds"]
            c = r.get("compute") or {}
            setup = c.get("setup_seconds", 0.0) if c else 0.0
            ax.bar(i, wall - setup, color=COLORS[b], width=0.7)
            if setup > 0:
                ax.bar(i, setup, bottom=wall - setup, color=COLORS[b], alpha=0.35, width=0.7, hatch="///",
                       edgecolor="white")
            txt = f"{wall:.1f} s" if wall < 100 else f"{wall:.0f} s"
            if base and b != "cpu":
                txt += f"\n{base / wall:.2g}×"
            ax.text(i, wall, txt, ha="center", va="bottom", fontsize=7.5, color=INK)
            xs.append(i)
            labels.append({"cpu": "CPU", "opencl:double": "OCL\ndbl", "opencl:mixed": "OCL\nmix",
                           "hip:double": "HIP\ndbl", "hip:mixed": "HIP\nmix"}[b])
        ax.set_xticks(xs)
        ax.set_xticklabels(labels, fontsize=7.5)
        ax.set_title(f"{CASE_LABELS.get(key[0], key[0])}\n{key[1]:,} ions", fontsize=9.5, color=INK)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.18)
        ax.grid(axis="x", visible=False)
    np.atleast_1d(axes)[0].set_ylabel("wall time of amcset run (s)")
    fig.suptitle("Wall time on a Ryzen 7 7800X3D (CPU, 16 threads) and an RX 7900 XTX (shared with other work: "
                 "indications). OCL = OpenCL, dbl/mix = double/mixed precision.\nHatched: setup (kernel compile, "
                 "upload); ×: speed-up over the CPU (below 1: slower)", color=INK, fontsize=10.5, y=1.07)
    fig.tight_layout()
    save(fig, out, "speed.svg")


def fig_occupancy(data, out):
    """Kernel steps per second against the lanes a run fills (one block of 16 ions per lane)."""
    p = os.path.join(data, "speed", "speed.json")
    if not os.path.exists(p):
        return
    rows = [r for r in json.load(open(p)) if r.get("compute") and r["compute"].get("kernel_seconds")]
    fig, ax = plt.subplots(figsize=(7.5, 4))
    for b in ["opencl:double", "opencl:mixed", "hip:double", "hip:mixed"]:
        pts = sorted((r["compute"]["lanes"], r["compute"]["total_steps"] / r["compute"]["kernel_seconds"], r["case"])
                     for r in rows if r["backend"] == b)
        if not pts:
            continue
        ax.plot([q[0] for q in pts], [q[1] for q in pts], "o", color=COLORS[b], label=LABELS[b], markersize=7,
                linestyle="none")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("lanes (work items) the run fills: blocks of 16 ions in flight")
    ax.set_ylabel("kernel steps per second")
    ax.axvspan(100, 6000, color="#fef2f2", zorder=0)
    ax.text(110, ax.get_ylim()[1] * 0.6, "too few blocks to fill the GPU:\nthe CPU is faster", fontsize=8.5,
            color="#b91c1c", va="top")
    ax.legend(loc="lower right")
    ax.set_title("One lane per block: a GPU needs about 10⁴ blocks in flight (RX 7900 XTX, every case of the "
                 "speed runs)", fontsize=10, color=INK, loc="left")
    save(fig, out, "occupancy.svg")


def fig_profile(out):
    """Numbers measured earlier and quoted in docs/GPU.md "Profile" (RX 7900 XTX, OpenCL, 50 keV Fe in Fe)."""
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.6), gridspec_kw={"width_ratios": [1.3, 1, 1]})
    ax = axes[0]
    labels = ["full kernel", "no tally\nwrites", "screened Rutherford\ninstead of magic", "all three\nsimplified"]
    vals = [6.1, 3.7, 2.5, 1.4]
    ax.barh(range(4), vals, color=["#ea580c", "#fb923c", "#fdba74", "#fed7aa"], height=0.6)
    for i, v in enumerate(vals):
        ax.text(v + 0.08, i, f"{v} ns", va="center", fontsize=9, color=INK)
    ax.set_yticks(range(4))
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("device time per kernel step (ns)")
    ax.set_title("Ablation (AMCSET_KERNEL_DEFINES)", fontsize=10, color=INK)
    ax.grid(axis="y", visible=False)
    ax = axes[1]
    ax.bar([0, 1], [5141, 188], color=["#2563eb", "#f59e0b"], width=0.6)
    for i, v in enumerate([5141, 188]):
        ax.text(i, v, f"{v}", ha="center", va="bottom", fontsize=9, color=INK)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["double-precision\ninstructions", "single-precision\ninstructions"], fontsize=8.5)
    ax.set_title("Radeon GPU Analyzer, gfx1100, mixed kernel\n(before moving the hot path to float)",
                 fontsize=10, color=INK)
    ax.grid(axis="x", visible=False)
    ax = axes[2]
    e = [10, 50, 200, 1000]
    d = [73, 196, 593, 1778]
    ax.plot(e, d, "o-", color="#dc2626")
    for x, y in zip(e, d):
        ax.text(x, y * 1.12, str(y), ha="center", fontsize=8.5, color=INK)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Fe ion energy in Fe (keV)")
    ax.set_ylabel("deepest recoil slot")
    ax.set_title("TRIM's recoil stack depth (400 ions each)", fontsize=10, color=INK)
    fig.tight_layout()
    save(fig, out, "profile.svg")


def anim_ring(out):
    """Batches (the first design) against the ring: lanes over time. An illustration."""
    rng = np.random.default_rng(3)
    lanes, blocks = 12, 72
    # Block lengths with the measured spread of full cascades (mean 118,000 steps, max 135,000), scaled
    # to 1, and a few long blocks as at 1 MeV.
    length = np.clip(rng.normal(1.0, 0.06, blocks), 0.85, 1.2)
    length[rng.choice(blocks, 4, replace=False)] *= 2.2
    batch = 24

    def batches():
        events = []  # (lane, block, start, end)
        t0 = 0.0
        for first in range(0, blocks, batch):
            ids = list(range(first, min(blocks, first + batch)))
            free = [t0] * lanes
            for b in ids:
                lane = int(np.argmin(free))
                events.append((lane, b, free[lane], free[lane] + length[b]))
                free[lane] += length[b]
            t0 = max(free)  # the batch ends with its slowest lane; then the host merges
            t0 += 0.15
        return events

    def ring(size=18):
        events = []
        free = [0.0] * lanes
        done = {}
        merged = 0
        for b in range(blocks):
            lane = int(np.argmin(free))
            start = free[lane]
            # A block needs a free slot: block b - size must be merged (all blocks before it finished).
            if b >= size:
                while merged <= b - size:
                    merged += 1
                start = max(start, max(done[j] for j in range(0, b - size + 1)))
            events.append((lane, b, start, start + length[b]))
            done[b] = start + length[b]
            free[lane] = start + length[b]
        return events

    ev = [batches(), ring()]
    t_end = max(e[3] for e in ev[0] + ev[1]) * 1.02
    fig, axes = plt.subplots(2, 1, figsize=(10, 5.4), sharex=True)
    titles = ["Batches (the first design): a batch lasts as long as its slowest lane",
              "Ring of block slots (now): a free lane takes the next block at once"]
    cmap = plt.get_cmap("tab20")
    frames = 90
    bars = []
    for ax, title, events in zip(axes, titles, ev):
        ax.set_xlim(0, t_end)
        ax.set_ylim(-0.6, lanes - 0.4)
        ax.set_yticks(range(lanes))
        ax.set_yticklabels([f"lane {i}" for i in range(lanes)], fontsize=7)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.grid(axis="y", visible=False)
        rects = []
        for lane, b, s, e in events:
            r = ax.barh(lane, 0, left=s, height=0.72, color=cmap(b % 20), edgecolor="white", linewidth=0.6)[0]
            rects.append((r, s, e, b))
        bars.append(rects)
    axes[1].set_xlabel("time (blocks of average length)")
    texts = [ax.text(1.0, 1.02, "", transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=INK)
             for ax in axes]
    lines = [ax.axvline(0, color=INK, linewidth=1) for ax in axes]

    def update(k):
        t = t_end * (k + 1) / frames
        for rects, txt, line, events in zip(bars, texts, lines, ev):
            busy = 0.0
            for r, s, e, b in rects:
                w = max(0.0, min(t, e) - s)
                r.set_width(w)
                busy += w
            finished = sum(1 for _, _, s, e in events if e <= t)
            util = busy / (lanes * t) if t > 0 else 0
            txt.set_text(f"blocks finished {finished}/{blocks}   lanes busy {100 * util:.0f} %")
            line.set_xdata([t, t])
        return [r for rects in bars for r, *_ in rects] + texts + lines

    fig.suptitle("Why the ring: the same 72 blocks on 12 lanes (illustration; block lengths with the measured "
                 "spread, four long ones)", color=INK, fontsize=10.5)
    fig.tight_layout()
    anim = animation.FuncAnimation(fig, update, frames=frames, interval=80, blit=False)
    p = os.path.join(out, "ring_vs_batches.gif")
    anim.save(p, writer=animation.PillowWriter(fps=12))
    plt.close(fig)
    print("wrote", p)


def anim_live(data, out):
    """What the live plots show during a GPU run: the first N ions (speed.py --live)."""
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "trim_compare"))
    import trimout  # noqa: E402

    dirs = sorted(glob.glob(os.path.join(data, "live", "[0-9]*")))
    frames = []
    for d in dirs:
        if not os.path.exists(os.path.join(d, "RANGE.txt")):
            continue
        r = trimout.read_range(os.path.join(d, "RANGE.txt"))
        v = trimout.read_vacancy(os.path.join(d, "VACANCY.txt"))
        s = json.load(open(os.path.join(d, "summary.json")))
        vac = np.sum(np.array(v["columns"]), axis=0) if v["columns"] else np.zeros(len(v["depth"]))
        frames.append({"n": int(s["ions"]), "depth": np.array(r["depth"]), "ions": np.array(r["ion_dist"]),
                       "vdepth": np.array(v["depth"]), "vac": vac, "range": s["range"]["mean_A"],
                       "strag": s["range"]["straggling_A"], "vpi": s["per_ion"]["vacancies"],
                       "device": (s.get("compute") or {}).get("device", "CPU")})
    if len(frames) < 2:
        return
    last = frames[-1]
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.9), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    a1, a2, a3 = axes
    norm = lambda f, key: f[key]  # noqa: E731  (TRIM's files are per ion already)
    a1.set_xlim(0, 600)
    a1.set_ylim(0, 1.15 * max(norm(f, "ions").max() for f in frames))
    a2.set_xlim(0, 600)
    a2.set_ylim(0, 1.15 * max(norm(f, "vac").max() for f in frames))
    (l1,) = a1.step(last["depth"], norm(last, "ions"), where="mid", color=COLORS["hip:mixed"])
    (l2,) = a2.step(last["vdepth"], norm(last, "vac"), where="mid", color="#dc2626")
    a1.set_title("Ion depth distribution", fontsize=10, color=INK)
    a1.set_xlabel("depth (Å)")
    a1.set_ylabel("(atoms/cm³) / (atoms/cm²)")
    a2.set_title("Vacancies (ion + recoils)", fontsize=10, color=INK)
    a2.set_xlabel("depth (Å)")
    a2.set_ylabel("vacancies / (Å · ion)")
    ns = np.array([f["n"] for f in frames])
    rng = np.array([f["range"] for f in frames])
    se = np.array([f["strag"] / math.sqrt(f["n"]) for f in frames])
    a3.set_xscale("log")
    a3.set_xlim(ns[0] * 0.8, ns[-1] * 1.25)
    a3.set_ylim(rng[-1] - 4 * se[0], rng[-1] + 4 * se[0])
    a3.axhline(rng[-1], color=MUTED, linewidth=0.8, linestyle=":")
    a3.set_title("Mean range ± standard error", fontsize=10, color=INK)
    a3.set_xlabel("ions merged")
    a3.set_ylabel("Å")
    eb = a3.errorbar([], [], yerr=[], fmt="o", color=COLORS["hip:mixed"], capsize=3)
    txt = fig.text(0.5, 0.94, "", ha="center", color=INK, fontsize=10.5)
    hold = 6

    def update(k):
        f = frames[min(k, len(frames) - 1)]
        l1.set_data(f["depth"], norm(f, "ions"))
        l2.set_data(f["vdepth"], norm(f, "vac"))
        upto = min(k, len(frames) - 1) + 1
        nonlocal eb
        eb.remove()
        eb = a3.errorbar(ns[:upto], rng[:upto], yerr=se[:upto], fmt="o", color=COLORS["hip:mixed"], capsize=3,
                         markersize=4)
        txt.set_text(f"50 keV Fe in Fe, full cascades, on {f['device'].split('[')[0].strip()}: "
                     f"{f['n']:,} ions merged — {f['vpi']:.1f} vacancies per ion")
        return l1, l2, txt

    fig.subplots_adjust(top=0.82, wspace=0.32, left=0.06, right=0.98, bottom=0.15)
    anim = animation.FuncAnimation(fig, update, frames=len(frames) + hold, interval=400, blit=False)
    p = os.path.join(out, "live_plots.gif")
    anim.save(p, writer=animation.PillowWriter(fps=3))
    plt.close(fig)
    print("wrote", p)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", default="build/gpu_compare")
    ap.add_argument("--out", default="docs/gpu/figures")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    style()
    fig_z_distribution(a.data, a.out)
    fig_z_by_case(a.data, a.out)
    fig_checks(a.data, a.out)
    fig_speed(a.data, a.out)
    fig_occupancy(a.data, a.out)
    fig_profile(a.out)
    anim_ring(a.out)
    anim_live(a.data, a.out)


if __name__ == "__main__":
    sys.exit(main())

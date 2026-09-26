#!/usr/bin/env python3
"""
animate_cascade.py

Animates a collision-cascade / bombardment simulation (SRIM-like CSV output)
in 3D, drawing particle tracks as they occur in time, colored by knock-on
generation, with a synchronized side panel of diagnostic time-series and a
running time indicator.

Usage
-----
    python animate_cascade.py bombardment_0.csv
    python animate_cascade.py bombardment_0.csv --duration 8
    python animate_cascade.py bombardment_0.csv --generations 0-2
    python animate_cascade.py bombardment_0.csv --generations 4
    python animate_cascade.py bombardment_0.csv --generations all
    python animate_cascade.py bombardment_0.csv --mkv out.mkv
    python animate_cascade.py bombardment_0.csv --mkv out.mkv --fps 30

Interactive mode (default) opens a matplotlib window you can rotate/zoom
while it plays. --mkv renders to a Matroska video instead (requires ffmpeg).

Data model
----------
The CSV has one row per simulation step, for *all* particles interleaved,
with no explicit particle ID. A new particle (knock-on) is recorded as a
row that shares the exact (x, y, z, time) of the row that "created" it,
but with generation += 1. The parent particle's track then later resumes
(appears again) at its next generation-matching row. This script
reconstructs individual contiguous particle tracks from that interleaved
stream using a stack keyed by generation depth, before animating them.
"""

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter

COLUMN_NAMES = [
    "energy", "x", "y", "z", "vx", "vy", "vz", "znum", "mass", "gen", "time",
]

GEN_CMAP = plt.get_cmap("plasma")
MAX_COLOR_GEN = 8  # generations beyond this reuse the last color


# --------------------------------------------------------------------------
# Data loading & track reconstruction
# --------------------------------------------------------------------------

def load_data(path):
    df = pd.read_csv(path)
    df.columns = COLUMN_NAMES
    df["gen"] = df["gen"].astype(int)
    df = df.reset_index(drop=True)
    return df


def reconstruct_tracks(df):
    """
    Rebuild contiguous per-particle tracks from the interleaved row stream.

    Returns a list of dicts, one per particle/track:
        {
          'gen': int,
          'rows': [row indices into df, in order],
          'parent': track_id or None,
        }
    """
    tracks = []
    stack = {}  # generation -> currently-open track id at that generation

    gens = df["gen"].values
    for i in range(len(df)):
        g = int(gens[i])
        if g not in stack:
            parent_tid = stack.get(g - 1, None)
            tid = len(tracks)
            tracks.append({"gen": g, "rows": [i], "parent": parent_tid})
            stack[g] = tid
        else:
            tid = stack[g]
            tracks[tid]["rows"].append(i)
        # Any track deeper than the current generation is now definitively
        # finished (a shallower row can only appear after all its deeper
        # descendants in the DFS-ordered stream have been fully emitted).
        for gg in [k for k in stack if k > g]:
            del stack[gg]

    return tracks


def build_track_arrays(df, tracks):
    """Attach numpy arrays of x, y, z, t, energy to each track for fast slicing."""
    x = df["x"].values
    y = df["y"].values
    z = df["z"].values
    t = df["time"].values
    e = df["energy"].values
    for tr in tracks:
        idx = np.array(tr["rows"])
        tr["idx"] = idx
        tr["t"] = t[idx]
        tr["x"] = x[idx]
        tr["y"] = y[idx]
        tr["z"] = z[idx]
        tr["e"] = e[idx]
        tr["t_start"] = tr["t"][0]
        tr["t_end"] = tr["t"][-1]
        tr["time_sorted"] = np.all(tr["t"][1:] >= tr["t"][:-1])
    return tracks


# --------------------------------------------------------------------------
# Derived diagnostic time-series
# --------------------------------------------------------------------------

def compute_diagnostics(df, tracks):
    """
    Compute, at every event (row) in time order:
      - cumulative displacement count (a knock-on event = 1 new displaced atom)
      - incident-particle (track 0) energy vs time
      - nuclear energy loss vs time: energy carried off by each newly created
        knock-on (gen > 0) at the moment it is created
      - electronic energy loss vs time: for *every* particle step (all
        generations), E_before - E_after - E_given_to_any_new_knockon
        produced at that step.

    Returns a dict of sorted (time, cumulative_value) arrays plus raw
    per-event arrays for reference.
    """
    t_all = df["time"].values
    e_all = df["energy"].values
    gen_all = df["gen"].values.astype(int)

    order = np.argsort(t_all, kind="mergesort")

    # ---- Displacement count: every row with generation > 0 represents a
    # newly displaced (knocked-on) atom created at that time.
    disp_times = t_all[gen_all > 0]
    disp_order = np.argsort(disp_times)
    disp_t_sorted = disp_times[disp_order]
    disp_cum = np.arange(1, len(disp_t_sorted) + 1)

    # ---- Incident particle (track index 0, generation 0) energy vs time
    incident = tracks[0]
    incident_t = incident["t"]
    incident_e = incident["e"]

    # ---- Nuclear energy loss vs time: energy of each new knock-on at its
    # creation time (== energy of first row of every track with gen>0)
    nuclear_t = []
    nuclear_e = []
    for tr in tracks:
        if tr["gen"] > 0:
            nuclear_t.append(tr["t"][0])
            nuclear_e.append(tr["e"][0])
    nuclear_t = np.array(nuclear_t)
    nuclear_e = np.array(nuclear_e)
    n_order = np.argsort(nuclear_t)
    nuclear_t_sorted = nuclear_t[n_order]
    nuclear_cum = np.cumsum(nuclear_e[n_order])

    # ---- Electronic energy loss vs time: for every step i -> i+1 within a
    # track (all generations), the energy is E_before at row i. Any
    # knock-ons spawned "at" a given row consume some energy that is
    # subtracted from what would otherwise go to nuclear stopping. Here we
    # follow the requested definition:
    #   electronic loss at a step = E_before - E_after - E_to_new_knockons
    # where E_before/E_after are consecutive energies of a single particle
    # along its own track, and E_to_new_knockons is the summed starting
    # energy of any knock-on(s) whose creation point/time coincides with
    # that step's end.
    #
    # Knock-ons are created at a specific (x,y,z,time) that matches the
    # *end* of some step of some track (the row immediately preceding the
    # knock-on's first row shares its position/time). We build a lookup
    # from (row index in df that "births" a knock-on) -> summed knock-on
    # energy, using the parent's row index recorded during reconstruction.
    birth_energy_by_parent_row = {}
    for tr in tracks:
        if tr["parent"] is not None:
            parent_tr = tracks[tr["parent"]]
            # the parent's row immediately preceding this child's first row,
            # i.e. the last parent row with the same time/position at the
            # branch point, is the row just before the child's row in the
            # dataframe ordering (== the birth row / trigger row).
            child_row0 = tr["rows"][0]
            trigger_row = child_row0 - 1
            birth_energy_by_parent_row[trigger_row] = (
                birth_energy_by_parent_row.get(trigger_row, 0.0) + tr["e"][0]
            )

    elec_t = []
    elec_loss = []
    for tr in tracks:
        rows = tr["rows"]
        es = tr["e"]
        ts = tr["t"]
        for k in range(len(rows) - 1):
            e_before = es[k]
            e_after = es[k + 1]
            df_row_of_before = rows[k]
            e_to_knockons = birth_energy_by_parent_row.get(df_row_of_before, 0.0)
            loss = e_before - e_after - e_to_knockons
            elec_t.append(ts[k + 1])
            elec_loss.append(max(loss, 0.0))
    elec_t = np.array(elec_t)
    elec_loss = np.array(elec_loss)
    e_order = np.argsort(elec_t)
    elec_t_sorted = elec_t[e_order]
    elec_cum = np.cumsum(elec_loss[e_order])

    return {
        "disp_t": disp_t_sorted, "disp_cum": disp_cum,
        "incident_t": incident_t, "incident_e": incident_e,
        "nuclear_t": nuclear_t_sorted, "nuclear_cum": nuclear_cum,
        "elec_t": elec_t_sorted, "elec_cum": elec_cum,
    }


# --------------------------------------------------------------------------
# Generation filtering
# --------------------------------------------------------------------------

def parse_generation_filter(spec, max_gen):
    """
    spec: 'all', a single int like '4', or a range like '0-2'.
    Returns a set of generation ints to display.
    """
    spec = spec.strip().lower()
    if spec == "all":
        return set(range(0, max_gen + 1))
    if "-" in spec:
        lo, hi = spec.split("-")
        lo, hi = int(lo), int(hi)
        return set(range(lo, hi + 1))
    return {int(spec)}


def gen_color(g):
    frac = min(g, MAX_COLOR_GEN) / MAX_COLOR_GEN
    return GEN_CMAP(0.15 + 0.85 * frac)


# --------------------------------------------------------------------------
# Animation
# --------------------------------------------------------------------------

def run(args):
    df = load_data(args.csv)
    tracks = reconstruct_tracks(df)
    tracks = build_track_arrays(df, tracks)
    diag = compute_diagnostics(df, tracks)

    max_gen = df["gen"].max()
    show_gens = parse_generation_filter(args.generations, max_gen)

    t_min = df["time"].min()
    data_t_max = df["time"].max()
    t_max = min(data_t_max, args.max_time) if args.max_time is not None else data_t_max
    if t_max <= t_min:
        raise ValueError(f"--max-time must be greater than the first event time ({t_min:.6g} s).")
    t_span = t_max - t_min

    # spatial bounds (small padding)
    pad = 0.05
    plotted_rows = df[df["time"] <= t_max]
    xr = (plotted_rows["x"].min(), plotted_rows["x"].max())
    yr = (plotted_rows["y"].min(), plotted_rows["y"].max())
    zr = (plotted_rows["z"].min(), plotted_rows["z"].max())

    def padded(r):
        span = r[1] - r[0]
        if span == 0:
            span = 1.0
        return (r[0] - pad * span, r[1] + pad * span)

    xr, yr, zr = padded(xr), padded(yr), padded(zr)

    # Only animate tracks whose generation is in the requested set
    visible_tracks = [tr for tr in tracks if tr["gen"] in show_gens]

    n_frames = (max(1, round(args.duration * args.fps))
                if args.mkv and args.fps is not None else args.frames)
    frame_times = np.linspace(t_min, t_max, n_frames)

    # ---- large 3D and XY views above full-width XZ and YZ views
    fig = plt.figure(figsize=(18, 12))
    fig.suptitle("Collision Cascade — Particle Tracks by Generation", fontsize=13, fontweight="bold")
    gs = GridSpec(6, 3, figure=fig, height_ratios=[1, 1, 1, 1, 1.2, 1.2],
                  width_ratios=[2.4, 1.05, 2], hspace=0.6, wspace=0.35,
                  left=0.04, right=0.98, top=0.90, bottom=0.08)

    ax3d = fig.add_subplot(gs[:4, 0], projection="3d")
    ax_xy = fig.add_subplot(gs[:4, 1])
    ax_xz = fig.add_subplot(gs[4, :])
    ax_yz = fig.add_subplot(gs[5, :])
    ax_disp = fig.add_subplot(gs[0, 2])
    ax_inc = fig.add_subplot(gs[1, 2])
    ax_nuc = fig.add_subplot(gs[2, 2])
    ax_ele = fig.add_subplot(gs[3, 2])

    ax3d.set_xlim(*xr)
    ax3d.set_ylim(*yr)
    ax3d.set_zlim(*zr)
    ax3d.set_xlabel("x (Å)")
    ax3d.set_ylabel("y (Å)")
    ax3d.set_zlabel("z (Å)")
    ax3d.set_title(f"Generations shown: {args.generations}", fontsize=10)

    ax_xy.set_xlim(*xr)
    ax_xy.set_ylim(*yr)
    ax_xy.set_aspect("equal", adjustable="box")
    ax_xy.set_xlabel("x (Å)")
    ax_xy.set_ylabel("y (Å)")
    ax_xy.set_title("XY view", fontsize=10)

    ax_xz.set_xlim(*zr)
    ax_xz.set_ylim(*xr)
    ax_xz.set_xlabel("z (Å)")
    ax_xz.set_ylabel("x (Å)")
    ax_xz.set_title("XZ view (Z horizontal)", fontsize=10)

    ax_yz.set_xlim(*zr)
    ax_yz.set_ylim(*yr)
    ax_yz.set_xlabel("z (Å)")
    ax_yz.set_ylabel("y (Å)")
    ax_yz.set_title("YZ view (Z horizontal)", fontsize=10)

    time_text = ax3d.text2D(0.02, 0.97, "", transform=ax3d.transAxes,
                             fontsize=11, va="top",
                             bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    # legend for generations present
    present_gens = sorted(show_gens & set(tr["gen"] for tr in visible_tracks))
    legend_handles = [Line2D([0], [0], color=gen_color(g), lw=2.5, label=f"Gen {g}")
                       for g in present_gens]
    if legend_handles:
        ax3d.legend(handles=legend_handles, loc="upper right", fontsize=8, framealpha=0.85)

    # line + head-marker artists per visible track
    line_artists = []
    head_artists = []
    xy_line_artists = []
    xy_head_artists = []
    xz_line_artists = []
    xz_head_artists = []
    yz_line_artists = []
    yz_head_artists = []
    for tr in visible_tracks:
        color = gen_color(tr["gen"])
        (line,) = ax3d.plot([], [], [], color=color, lw=1.2, alpha=0.85)
        (head,) = ax3d.plot([], [], [], marker="o", color=color, markersize=3.5,
                             markeredgecolor="black", markeredgewidth=0.3)
        line_artists.append(line)
        head_artists.append(head)
        (xy_line,) = ax_xy.plot([], [], color=color, lw=1.2, alpha=0.85)
        (xy_head,) = ax_xy.plot([], [], marker="o", color=color, markersize=3.5,
                                markeredgecolor="black", markeredgewidth=0.3)
        xy_line_artists.append(xy_line)
        xy_head_artists.append(xy_head)
        (xz_line,) = ax_xz.plot([], [], color=color, lw=1.2, alpha=0.85)
        (xz_head,) = ax_xz.plot([], [], marker="o", color=color, markersize=3.5,
                                markeredgecolor="black", markeredgewidth=0.3)
        xz_line_artists.append(xz_line)
        xz_head_artists.append(xz_head)
        (yz_line,) = ax_yz.plot([], [], color=color, lw=1.2, alpha=0.85)
        (yz_head,) = ax_yz.plot([], [], marker="o", color=color, markersize=3.5,
                    markeredgecolor="black", markeredgewidth=0.3)
        yz_line_artists.append(yz_line)
        yz_head_artists.append(yz_head)

    # ---- diagnostic panels: static full curves, with a moving cursor
    def setup_panel(ax, xdata, ydata, ylabel, color, scientific=False):
        visible = xdata <= t_max
        xdata = xdata[visible]
        ydata = ydata[visible]
        ax.plot(xdata, ydata, color=color, lw=1.3)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(labelsize=7)
        if scientific:
            formatter = ScalarFormatter(useMathText=True)
            formatter.set_scientific(True)
            formatter.set_powerlimits((0, 0))
            ax.yaxis.set_major_formatter(formatter)
        ax.set_xlim(t_min, t_max)
        ymax = ydata.max() if len(ydata) else 1.0
        ax.set_ylim(0, ymax * 1.1 if ymax > 0 else 1.0)
        (cursor,) = ax.plot([t_min], [ydata[0] if len(ydata) else 0], "o", color="red",
                             markersize=5, zorder=5)
        vline = ax.axvline(t_min, color="red", lw=0.8, alpha=0.6)
        return cursor, vline

    cur_disp, vl_disp = setup_panel(ax_disp, diag["disp_t"], diag["disp_cum"],
                                     "# Displacements", "#2166ac")
    ax_disp.set_title("Cumulative displacements vs time", fontsize=9)

    cur_inc, vl_inc = setup_panel(ax_inc, diag["incident_t"], diag["incident_e"],
                                   "Incident E (eV)", "#b2182b", scientific=True)
    ax_inc.set_title("Incident particle energy vs time", fontsize=9)

    cur_nuc, vl_nuc = setup_panel(ax_nuc, diag["nuclear_t"], diag["nuclear_cum"],
                                   "Nuclear loss (eV)", "#4daf4a", scientific=True)
    ax_nuc.set_title("Cumulative nuclear energy loss vs time", fontsize=9)

    cur_ele, vl_ele = setup_panel(ax_ele, diag["elec_t"], diag["elec_cum"],
                                   "Electronic loss (eV)", "#984ea3", scientific=True)
    ax_ele.set_title("Cumulative electronic energy loss vs time", fontsize=9)
    ax_ele.set_xlabel("time (s)", fontsize=8)

    def nearest_cum_value(xs, ys, tcur):
        """Value of a cumulative step function ys at time tcur (last point <= tcur)."""
        if len(xs) == 0:
            return 0.0
        pos = np.searchsorted(xs, tcur, side="right") - 1
        if pos < 0:
            return 0.0
        return ys[min(pos, len(ys) - 1)]

    def update_panel(cursor, vline, xs, ys, tcur):
        val = nearest_cum_value(xs, ys, tcur)
        cursor.set_data([tcur], [val])
        vline.set_xdata([tcur, tcur])

    def fmt_time(t):
        return f"t = {t:.3e} s"

    def init():
        for line, head in zip(line_artists, head_artists):
            line.set_data([], [])
            line.set_3d_properties([])
            head.set_data([], [])
            head.set_3d_properties([])
        for line, head in zip(xy_line_artists, xy_head_artists):
            line.set_data([], [])
            head.set_data([], [])
        for line, head in zip(xz_line_artists, xz_head_artists):
            line.set_data([], [])
            head.set_data([], [])
        for line, head in zip(yz_line_artists, yz_head_artists):
            line.set_data([], [])
            head.set_data([], [])
        time_text.set_text(fmt_time(t_min))
        return (line_artists + head_artists + xy_line_artists + xy_head_artists
                + xz_line_artists + xz_head_artists + yz_line_artists
                + yz_head_artists + [time_text])

    def animate(frame_i):
        tcur = frame_times[frame_i]
        for tr, line, head, xy_line, xy_head, xz_line, xz_head, yz_line, yz_head in zip(
                visible_tracks, line_artists, head_artists, xy_line_artists,
            xy_head_artists, xz_line_artists, xz_head_artists,
            yz_line_artists, yz_head_artists):
            if tcur < tr["t_start"]:
                line.set_data([], [])
                line.set_3d_properties([])
                head.set_data([], [])
                head.set_3d_properties([])
                xy_line.set_data([], [])
                xy_head.set_data([], [])
                xz_line.set_data([], [])
                xz_head.set_data([], [])
                yz_line.set_data([], [])
                yz_head.set_data([], [])
                continue
            if tr["time_sorted"]:
                n = np.searchsorted(tr["t"], tcur, side="right")
                if n == 0:
                    continue
                xs = tr["x"][:n]
                ys = tr["y"][:n]
                zs = tr["z"][:n]
            else:
                visible = tr["t"] <= tcur
                if not visible.any():
                    continue
                xs = tr["x"][visible]
                ys = tr["y"][visible]
                zs = tr["z"][visible]
            if len(xs) == 0:
                continue
            line.set_data(xs, ys)
            line.set_3d_properties(zs)
            head.set_data([xs[-1]], [ys[-1]])
            head.set_3d_properties([zs[-1]])
            xy_line.set_data(xs, ys)
            xy_head.set_data([xs[-1]], [ys[-1]])
            xz_line.set_data(zs, xs)
            xz_head.set_data([zs[-1]], [xs[-1]])
            yz_line.set_data(zs, ys)
            yz_head.set_data([zs[-1]], [ys[-1]])

        time_text.set_text(fmt_time(tcur))

        update_panel(cur_disp, vl_disp, diag["disp_t"], diag["disp_cum"], tcur)
        update_panel(cur_inc, vl_inc, diag["incident_t"], diag["incident_e"], tcur)
        update_panel(cur_nuc, vl_nuc, diag["nuclear_t"], diag["nuclear_cum"], tcur)
        update_panel(cur_ele, vl_ele, diag["elec_t"], diag["elec_cum"], tcur)

        return (line_artists + head_artists + xy_line_artists + xy_head_artists
            + xz_line_artists + xz_head_artists + yz_line_artists + yz_head_artists
            + [time_text, cur_disp, cur_inc, cur_nuc, cur_ele])

    interval_ms = (args.duration * 1000.0) / n_frames

    ani = animation.FuncAnimation(
        fig, animate, frames=n_frames, init_func=init,
        interval=interval_ms, blit=False, repeat=True,
    )

    if args.mkv:
        fps = args.fps if args.fps else max(1, round(n_frames / args.duration))
        writer = animation.FFMpegWriter(fps=fps, bitrate=args.bitrate,
                                        extra_args=["-f", "matroska"])
        playback_duration = n_frames / fps
        print(f"Rendering MKV to {args.mkv} ({n_frames} frames @ {fps} fps, "
              f"~{playback_duration:.1f}s playback)...")
        ani.save(args.mkv, writer=writer, dpi=args.dpi)
        print("Done.")
    else:
        plt.show()


def build_arg_parser():
    p = argparse.ArgumentParser(description="Animate a bombardment/cascade CSV in 3D.")
    p.add_argument("csv", help="Path to bombardment CSV file.")
    p.add_argument("--duration", type=float, default=5.0,
                    help="Total animation playback duration in seconds (default: 5).")
    p.add_argument("--max-time", type=float, default=None,
                    help="Latest simulation time to plot in seconds (default: all available data).")
    p.add_argument("--frames", type=int, default=300,
                    help="Frame count for interactive playback and video export when --fps is omitted (default: 300).")
    p.add_argument("--generations", type=str, default="all",
                    help="Which generations to show: 'all', a single number "
                         "(e.g. '2'), or a range (e.g. '0-3'). Default: all.")
    p.add_argument("--mkv", type=str, default=None,
                    help="If given, render to this Matroska file path instead of "
                         "opening an interactive window (requires ffmpeg).")
    p.add_argument("--fps", type=int, default=None,
                    help="Video frames per second; with --duration, determines export frame count and overrides --frames.")
    p.add_argument("--dpi", type=int, default=120, help="DPI for video export.")
    p.add_argument("--bitrate", type=int, default=3000, help="Bitrate (kbps) for video export.")
    return p


def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
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

# animate_cascade.jl
#
#   julia -t auto --project animate_cascade.jl bombardment_0.csv
#   julia -t auto --project animate_cascade.jl bombardment_0.csv --generations 0-2
#   julia -t auto --project animate_cascade.jl bombardment_0.csv --mkv out.mkv --fps 30

using ArgParse, CSV, GLMakie, Printf

const COLNAMES = [:energy, :x, :y, :z, :vx, :vy, :vz, :znum, :mass, :gen, :time]
const MAX_COLOR_GEN = 8
const GEN_CMAP = cgrad(:plasma)

gen_color(g) = RGBAf(GEN_CMAP[0.15 + 0.85 * min(g, MAX_COLOR_GEN) / MAX_COLOR_GEN])

# --------------------------------------------------------------------------
# Loading & track reconstruction
# --------------------------------------------------------------------------

function load_data(path)
    f = CSV.File(path; header=COLNAMES, skipto=2, types=Float64,
                 ntasks=Threads.nthreads())
    return (energy = Vector{Float64}(f.energy),
            x = Vector{Float64}(f.x), y = Vector{Float64}(f.y),
            z = Vector{Float64}(f.z), t = Vector{Float64}(f.time),
            gen = round.(Int, f.gen))
end

struct Tracks
    gen::Vector{Int}
    parent::Vector{Int32}
    start::Vector{Int}        # rows of track j: rows[start[j]:start[j+1]-1]
    rows::Vector{Int32}       # row indices grouped by track, in stream order
end

function reconstruct_tracks(gen::Vector{Int})
    n = length(gen)
    stack = zeros(Int32, maximum(gen) + 2)          # index g+1 -> open track id
    row_track = Vector{Int32}(undef, n)
    tgen = Int[]
    tparent = Int32[]
    @inbounds for i in 1:n
        g = gen[i]
        tid = stack[g + 1]
        if tid == 0
            push!(tgen, g)
            push!(tparent, g > 0 ? stack[g] : Int32(0))
            tid = Int32(length(tgen))
            stack[g + 1] = tid
        end
        row_track[i] = tid
        for k in (g + 2):length(stack)
            stack[k] = 0
        end
    end
    nt = length(tgen)
    cnt = zeros(Int, nt)
    @inbounds for tid in row_track
        cnt[tid] += 1
    end
    start = Vector{Int}(undef, nt + 1)
    start[1] = 1
    @inbounds for j in 1:nt
        start[j + 1] = start[j] + cnt[j]
    end
    ptr = start[1:nt]
    rows = Vector{Int32}(undef, n)
    @inbounds for i in 1:n
        tid = row_track[i]
        rows[ptr[tid]] = i
        ptr[tid] += 1
    end
    return Tracks(tgen, tparent, start, rows)
end

# --------------------------------------------------------------------------
# Diagnostics
# --------------------------------------------------------------------------

function compute_diagnostics(d, T::Tracks)
    t, e, gen = d.t, d.energy, d.gen
    nt = length(T.gen)

    disp_t = sort(t[gen .> 0])
    disp_cum = collect(1.0:length(disp_t))

    r1 = T.rows[T.start[1]:T.start[2]-1]
    incident_t, incident_e = t[r1], e[r1]

    nt_t = Float64[]; nt_e = Float64[]
    birth = zeros(length(t))
    for j in 1:nt
        r0 = T.rows[T.start[j]]
        if T.gen[j] > 0
            push!(nt_t, t[r0]); push!(nt_e, e[r0])
        end
        if T.parent[j] != 0
            birth[r0 - 1] += e[r0]
        end
    end
    p = sortperm(nt_t)
    nuclear_t = nt_t[p]
    nuclear_cum = cumsum(nt_e[p])

    elec_t = Vector{Float64}(undef, length(t) - nt)
    elec_l = similar(elec_t)
    m = 0
    @inbounds for j in 1:nt
        s, en = T.start[j], T.start[j + 1] - 1
        for k in s:en-1
            r, rn = T.rows[k], T.rows[k + 1]
            m += 1
            elec_t[m] = t[rn]
            elec_l[m] = max(e[r] - e[rn] - birth[r], 0.0)
        end
    end
    resize!(elec_t, m); resize!(elec_l, m)
    p = sortperm(elec_t)
    return (; disp_t, disp_cum, incident_t, incident_e,
            nuclear_t, nuclear_cum, elec_t = elec_t[p], elec_cum = cumsum(elec_l[p]))
end

function parse_generation_filter(spec, maxgen)
    spec = lowercase(strip(spec))
    mask = falses(maxgen + 1)
    if spec == "all"
        mask .= true
    elseif occursin("-", spec)
        lo, hi = parse.(Int, split(spec, "-"))
        for g in lo:min(hi, maxgen); mask[g + 1] = true; end
    else
        g = parse(Int, spec)
        g <= maxgen && (mask[g + 1] = true)
    end
    return mask
end

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

function run(args)
    d = load_data(args["csv"])
    T = reconstruct_tracks(d.gen)
    diag = compute_diagnostics(d, T)
    n = length(d.t)
    nt = length(T.gen)
    maxgen = maximum(d.gen)
    showmask = parse_generation_filter(args["generations"], maxgen)

    t_min = minimum(d.t)
    t_max = args["max-time"] === nothing ? maximum(d.t) : min(maximum(d.t), args["max-time"])
    t_max > t_min || error("--max-time must be greater than the first event time ($t_min s).")

    # spatial bounds
    sel = d.t .<= t_max
    function padded(v)
        lo, hi = extrema(view(v, sel))
        span = hi == lo ? 1.0 : hi - lo
        return (lo - 0.05span, hi + 0.05span)
    end
    xr, yr, zr = padded(d.x), padded(d.y), padded(d.z)

    # ---- flat vertex buffer: one polyline for ALL visible tracks, NaN-separated
    bufidx  = zeros(Int32, n)     # row -> vertex slot (0 = hidden)
    rowhead = zeros(Int32, n)     # row -> head-marker slot
    colors  = RGBAf[]
    hcolors = RGBAf[]
    pos = 0
    for j in 1:nt
        showmask[T.gen[j] + 1] || continue
        c = gen_color(T.gen[j])
        push!(hcolors, c)
        h = length(hcolors)
        for k in T.start[j]:T.start[j + 1]-1
            r = T.rows[k]
            pos += 1
            bufidx[r] = pos
            rowhead[r] = h
            push!(colors, c)
        end
        pos += 1; push!(colors, c)          # NaN separator slot
    end
    nbuf, nheads = pos, length(hcolors)

    # ---- global time-ordered event list (only rows that are actually drawn)
    order = sortperm(d.t; alg=MergeSort)
    ev = Int32[r for r in order if bufidx[r] != 0 && d.t[r] <= t_max]
    ev_t = d.t[ev]

    x32, y32, z32 = Float32.(d.x), Float32.(d.y), Float32.(d.z)
    nan3, nan2 = Point3f(NaN, NaN, NaN), Point2f(NaN, NaN)
    p3  = fill(nan3, nbuf);  pxy = fill(nan2, nbuf)
    pxz = fill(nan2, nbuf);  pyz = fill(nan2, nbuf)
    h3  = fill(nan3, nheads); hxy = fill(nan2, nheads)
    hxz = fill(nan2, nheads); hyz = fill(nan2, nheads)
    o3, oxy, oxz, oyz = Observable(p3), Observable(pxy), Observable(pxz), Observable(pyz)
    oh3, ohxy, ohxz, ohyz = Observable(h3), Observable(hxy), Observable(hxz), Observable(hyz)

    # ---- figure
    fig = Figure(size=(1800, 1200))
    title_obs = Observable(@sprintf("Collision Cascade — Particle Tracks by Generation     t = %.3e s", t_min))
    Label(fig[0, 1:3], title_obs; fontsize=18, font=:bold)

    ax3d = Axis3(fig[1, 1]; xlabel="x (Å)", ylabel="y (Å)", zlabel="z (Å)",
                 title="Generations shown: $(args["generations"])", aspect=(4, 4, 3))
    ax_xy = Axis(fig[1, 2]; xlabel="x (Å)", ylabel="y (Å)", title="XY view", aspect=DataAspect())
    ax_xz = Axis(fig[2, 1:3]; xlabel="z (Å)", ylabel="x (Å)", title="XZ view (Z horizontal)")
    ax_yz = Axis(fig[3, 1:3]; xlabel="z (Å)", ylabel="y (Å)", title="YZ view (Z horizontal)")
    limits!(ax3d, xr..., yr..., zr...)
    limits!(ax_xy, xr..., yr...)
    limits!(ax_xz, zr..., xr...)
    limits!(ax_yz, zr..., yr...)

    for (ax, o, oh) in ((ax3d, o3, oh3), (ax_xy, oxy, ohxy), (ax_xz, oxz, ohxz), (ax_yz, oyz, ohyz))
        lines!(ax, o; color=colors, linewidth=1.2)
        scatter!(ax, oh; color=hcolors, markersize=6, strokewidth=0.5, strokecolor=:black)
    end

    present = sort(unique(T.gen[[showmask[g + 1] for g in T.gen]]))
    isempty(present) || axislegend(ax3d, [LineElement(color=gen_color(g), linewidth=3) for g in present],
                                   ["Gen $g" for g in present]; position=:rt, labelsize=11)

    # ---- diagnostic panels
    gd = fig[1, 3] = GridLayout()
    sci(vals) = [@sprintf("%.1e", v) for v in vals]
    function setup_panel(row, xs, ys, ylabel, ttl, color; scientific=false, xlabel="")
        ax = Axis(gd[row, 1]; ylabel, title=ttl, xlabel, titlesize=12, ylabelsize=11,
                  ytickformat = scientific ? sci : Makie.automatic)
        keep = xs .<= t_max
        xs, ys = xs[keep], ys[keep]
        lines!(ax, xs, ys; color, linewidth=1.3)
        ymax = isempty(ys) ? 1.0 : maximum(ys)
        xlims!(ax, t_min, t_max); ylims!(ax, 0, ymax > 0 ? 1.1ymax : 1.0)
        cur = Observable(Point2f(t_min, isempty(ys) ? 0f0 : ys[1]))
        vl = Observable(t_min)
        vlines!(ax, vl; color=(:red, 0.6), linewidth=1)
        scatter!(ax, cur; color=:red, markersize=8)
        return (cur, vl, xs, ys)
    end
    panels = [
        setup_panel(1, diag.disp_t, diag.disp_cum, "# Displacements", "Cumulative displacements vs time", "#2166ac"),
        setup_panel(2, diag.incident_t, diag.incident_e, "Incident E (eV)", "Incident particle energy vs time", "#b2182b"; scientific=true),
        setup_panel(3, diag.nuclear_t, diag.nuclear_cum, "Nuclear loss (eV)", "Cumulative nuclear energy loss vs time", "#4daf4a"; scientific=true),
        setup_panel(4, diag.elec_t, diag.elec_cum, "Electronic loss (eV)", "Cumulative electronic energy loss vs time", "#984ea3"; scientific=true, xlabel="time (s)"),
    ]

    colsize!(fig.layout, 1, Auto(2.4)); colsize!(fig.layout, 2, Auto(1.05)); colsize!(fig.layout, 3, Auto(2.0))
    rowsize!(fig.layout, 1, Auto(4.0)); rowsize!(fig.layout, 2, Auto(1.2)); rowsize!(fig.layout, 3, Auto(1.2))

    # ---- incremental time stepping: only newly revealed events are written
    applied = Ref(0)
    function reset!()
        fill!(p3, nan3); fill!(pxy, nan2); fill!(pxz, nan2); fill!(pyz, nan2)
        fill!(h3, nan3); fill!(hxy, nan2); fill!(hxz, nan2); fill!(hyz, nan2)
        applied[] = 0
    end
    cumval(xs, ys, t) = (p = searchsortedlast(xs, t); p < 1 ? 0f0 : Float32(ys[p]))

    function advance_to!(tcur)
        target = searchsortedlast(ev_t, tcur)
        target < applied[] && reset!()
        if target > applied[]
            @inbounds for k in applied[]+1:target
                r = ev[k]; b = bufidx[r]; h = rowhead[r]
                xx, yy, zz = x32[r], y32[r], z32[r]
                p3[b] = Point3f(xx, yy, zz); pxy[b] = Point2f(xx, yy)
                pxz[b] = Point2f(zz, xx);    pyz[b] = Point2f(zz, yy)
                h3[h] = p3[b]; hxy[h] = pxy[b]; hxz[h] = pxz[b]; hyz[h] = pyz[b]
            end
            applied[] = target
            notify(o3); notify(oxy); notify(oxz); notify(oyz)
            notify(oh3); notify(ohxy); notify(ohxz); notify(ohyz)
        end
        title_obs[] = @sprintf("Collision Cascade — Particle Tracks by Generation     t = %.3e s", tcur)
        for (cur, vl, xs, ys) in panels
            cur[] = Point2f(tcur, cumval(xs, ys, tcur))
            vl[] = tcur
        end
        return nothing
    end

    duration = args["duration"]
    if args["mkv"] !== nothing
        fps = args["fps"]
        nframes = fps === nothing ? args["frames"] : max(1, round(Int, duration * fps))
        fps = fps === nothing ? max(1, round(Int, nframes / duration)) : fps
        frame_times = range(t_min, t_max; length=nframes)
        println("Rendering MKV to $(args["mkv"]) ($nframes frames @ $fps fps, ~$(round(nframes / fps; digits=1))s playback)...")
        record(fig, args["mkv"], 1:nframes; framerate=fps,
               px_per_unit=args["dpi"] / 100, compression=args["crf"]) do i
            advance_to!(frame_times[i])
        end
        println("Done.")
    else
        # ---- playback controls
        ctl = fig[4, 1:3] = GridLayout()
        playbtn = Button(ctl[1:3, 1]; label = "Pause", width = 80)
        Label(ctl[1, 2], "time";       halign = :right)
        Label(ctl[2, 2], "loop range"; halign = :right)
        scrub = Slider(ctl[1, 3]; range = 0:0.0005:1, startvalue = 0)
        rng   = IntervalSlider(ctl[2, 3]; range = 0:0.0005:1, startvalues = (0.0, 1.0))
        # speed slider works in log10 space: -2 -> 0.01x, 0 -> 1x, 2 -> 100x
        spd   = Slider(ctl[3, 3]; range = -2:0.02:2, startvalue = 0)
        Label(ctl[3, 2], lift(v -> @sprintf("speed %.2f×", 10.0^v), spd.value);
              halign = :right)
        colsize!(ctl, 3, Auto(1.0))

        playing  = Ref(true)
        pos      = Ref(0.0)          # playback position, fraction of [t_min, t_max]
        updating = Ref(false)        # guards against loop <-> slider feedback
        tof(f) = t_min + f * (t_max - t_min)

        toggle!() = (playing[] = !playing[];
                     playbtn.label[] = playing[] ? "Pause" : "Play")
        on(_ -> toggle!(), playbtn.clicks)
        on(events(fig).keyboardbutton) do ev
            ev.action == Keyboard.press || return
            if ev.key == Keyboard.space
                toggle!()
            elseif ev.key == Keyboard.up
                set_close_to!(spd, min(spd.value[] + log10(2), 2.0))
            elseif ev.key == Keyboard.down
                set_close_to!(spd, max(spd.value[] - log10(2), -2.0))
            end
        end
        on(scrub.value) do f            # user dragging the time slider
            updating[] && return
            pos[] = f
            advance_to!(tof(f))
        end

        screen = display(fig)
        last = time()
        while isopen(screen)
            now = time(); dt = now - last; last = now
            if playing[]
                lo, hi = rng.interval[]
                speed = 10.0^spd.value[]
                p = pos[] + dt * speed / duration
                (p >= hi || p < lo) && (p = lo)      # loop within the chosen range
                pos[] = p
                advance_to!(tof(p))
                updating[] = true
                set_close_to!(scrub, p)
                updating[] = false
            end
            sleep(1 / 120)
        end
    end
    return nothing
end

function main()
    s = ArgParseSettings(description="Animate a bombardment/cascade CSV in 3D.")
    @add_arg_table! s begin
        "csv"
            help = "Path to bombardment CSV file."
            required = true
        "--duration"
            arg_type = Float64
            default = 5.0
            help = "Total playback duration in seconds."
        "--max-time"
            arg_type = Float64
            default = nothing
            help = "Latest simulation time to plot in seconds."
        "--frames"
            arg_type = Int
            default = 300
            help = "Video frame count when --fps is omitted."
        "--generations"
            default = "all"
            help = "'all', a number (e.g. '2'), or a range (e.g. '0-3')."
        "--mkv"
            default = nothing
            help = "Render to this Matroska file instead of opening a window."
        "--fps"
            arg_type = Int
            default = nothing
        "--dpi"
            arg_type = Int
            default = 120
        "--crf"
            arg_type = Int
            default = 20
            help = "x264 quality (lower = better); replaces --bitrate."
    end
    run(parse_args(s))
end

main()
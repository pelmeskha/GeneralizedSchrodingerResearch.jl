module Plotting

using Plots

export plot_coupled_gl, plot_profiles, plot_error_profiles, plot_true_errors, plot_convergence, fit_order, plot_energy_exchange, plot_energy_evolution, plot_energy_total, envelope_clip_recorder, plot_gamma_grid


function save_plot(p, prefix, suffix, outdir)
    mkpath(outdir)
    path = joinpath(outdir, string(prefix, suffix))
    savefig(p, path)
    return path
end

"""
    plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis, outdir="results")

    Save two heatmaps showing the evolution of |u| and |v| as functions of z and T.

    - `T` : vector of retarded time values (length N_T)
    - `z` : vector of propagation distances (length N_z)
    - `U` : complex array of size (N_T, N_z) containing u(z,T)
    - `V` : complex array of size (N_T, N_z) containing v(z,T)
    - `prefix` : string prepended to output filenames
    - `outdir` : directory for saved figures (created if missing)
    - `cmap`  : colormap (symbol), e.g., :plasma, :inferno, :turbo
    - `max_z_points` : at most this many z-columns are drawn (every n-th column is taken);
      the saved figure cannot show more columns than its pixel width
"""
function plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis, outdir="results", max_z_points=2000)
    cols = 1:cld(length(z), max_z_points):length(z)
    z_plot = z[cols]

    # |u|
    p1 = heatmap(
        z_plot,
        T,
        abs.(U[:, cols]),
        xlabel = "z", ylabel = "t",
        color  = cmap,
        dpi    = 300
    )
    save_plot(p1, prefix, "_abs_u.png", outdir)

    # |v|
    p2 = heatmap(
        z_plot,
        T,
        abs.(V[:, cols]),
        xlabel = "z", ylabel = "t",
        color  = cmap,
        dpi    = 300
    )
    save_plot(p2, prefix, "_abs_v.png", outdir)
end
"""
    plot_coupled_gl(T, z, true_u::Function, true_v::Function; kwargs...)

    Evaluate analytical amplitudes on the (T, z) grid and plot their heatmaps.
    Analytical functions must have signature field(z, T).
"""
function plot_coupled_gl(T, z, true_u::Function, true_v::Function; kwargs...)
    U_true_abs = [abs(true_u(zi, tj)) for tj in T, zi in z]
    V_true_abs = [abs(true_v(zi, tj)) for tj in T, zi in z]

    return plot_coupled_gl(T, z, U_true_abs, V_true_abs; kwargs...)
end
const PROFILE_COLORS = [:blue, :red, :green, :orange, :purple, :brown]

format_z(z) = isinteger(z) ? string(Int(z)) : string(z)

function sample_profile(solution::Function, T_grid, z_grid, z_value)
    return [solution(z_value, t) for t in T_grid]
end

function sample_profile(solution::AbstractMatrix, T_grid, z_grid, z_value)
    size(solution) == (length(T_grid), length(z_grid)) ||
        throw(DimensionMismatch("solution must have size $((length(T_grid), length(z_grid)))"))
    j = argmin(abs.(z_grid .- z_value))
    abs(z_grid[j] - z_value) <= 1e-8 * max(1.0, abs(z_value)) ||
        throw(ArgumentError("z=$(format_z(z_value)) is not on the z grid"))
    return solution[:, j]
end

"""
    plot_profiles(; T_grid, z_grid, z_timestamps, u_solutions=nothing, v_solutions=nothing,
                  prefixes=nothing, prefix="profiles", dpi=300, outdir="results")

    Plot |u(T)| and |v(T)| profiles at every z in `z_timestamps`, one figure per z and component.
    Every solution is drawn on the same figure in its own color.

    - `T_grid` : retarded time array (length N_T)
    - `z_grid` : propagation distances of the numerical grid (length N_z)
    - `z_timestamps` : z values at which profiles are plotted; each must lie on `z_grid`
    - `u_solutions` : solutions for u. Each element is a complex matrix (N_T, N_z)
      or a callable `field(z, T)`
    - `v_solutions` : solutions for v in the same format. `nothing` by default, so no v plots.
      If both `u_solutions` and `v_solutions` are `nothing`, nothing is plotted
    - `prefixes` : one name per solution, used in legends as "u, <name>, at z=..."
      (defaults to "solution 1", "solution 2", ...)
    - `colors`, `linewidths`, `linestyles` : one value per solution. Default: the palette
      `PROFILE_COLORS`, width 1.5 for the first solution and 3.5 for the others, `:solid` for the
      first and `:dot` for the others
    - `legend` : `false` removes the legend
    - `prefix` : filename prefix for saved figures
    - `dpi` : resolution
    - `outdir` : directory for saved figures (created if missing)

    Returns the vector of plots.
"""
function plot_profiles(;
    T_grid, z_grid, z_timestamps,
    u_solutions=nothing, v_solutions=nothing,
    prefixes=nothing,
    colors=nothing, linewidths=nothing, linestyles=nothing, legend=true,
    prefix="profiles", dpi=300, outdir="results"
)
    u_list = isnothing(u_solutions) ? [] : u_solutions
    v_list = isnothing(v_solutions) ? [] : v_solutions
    n_solutions = max(length(u_list), length(v_list))
    names = isnothing(prefixes) ? ["solution $i" for i in 1:n_solutions] : prefixes
    length(names) >= n_solutions ||
        throw(ArgumentError("prefixes must contain at least $n_solutions names"))

    plots = []
    for (component, solutions) in (("u", u_list), ("v", v_list))
        isempty(solutions) && continue
        for z_value in z_timestamps
            z_label = format_z(z_value)
            p = plot(
                xlabel = "t",
                ylabel = "|$component|",
                legend = legend ? :best : false,
                dpi    = dpi
            )
            for (i, solution) in enumerate(solutions)
                profile = abs.(sample_profile(solution, T_grid, z_grid, z_value))
                plot!(
                    p, T_grid, profile;
                    label = "$component, $(names[i]), at z=$z_label",
                    linewidth = isnothing(linewidths) ? (i == 1 ? 1.5 : 3.5) : linewidths[i],
                    linestyle = isnothing(linestyles) ? (i == 1 ? :solid : :dot) : linestyles[i],
                    color = isnothing(colors) ? PROFILE_COLORS[mod1(i, length(PROFILE_COLORS))] : colors[i]
                )
            end
            save_plot(p, prefix, "_$(component)_z=$z_label.png", outdir)
            push!(plots, p)
        end
    end

    return plots
end
"""
    plot_error_profiles(T, z, U, V; prefix="error", dpi=300, outdir="results")

    Compute and plot the relative change of |u| and |v| from start to end.

    The relative error is defined as
        ( |field(T, z_end)| - |field(T, z_start)| ) / max_T |field(T, z_start)|
    so that the curve shows the change normalised to the initial peak amplitude.
"""
function plot_error_profiles(T, z, U, V; prefix="error", dpi=300, outdir="results")
    # Initial and final columns
    u_init = abs.(U[:, 1])
    u_final = abs.(U[:, end])
    v_init = abs.(V[:, 1])
    v_final = abs.(V[:, end])

    # Normalization constants (max of initial amplitude)
    u_scale = maximum(u_init)
    v_scale = maximum(v_init)

    # Relative error (guard against zero max, though never for physical pulses)
    u_err = u_scale == 0 ? zero(u_init) : (u_final .- u_init) ./ u_scale
    v_err = v_scale == 0 ? zero(v_init) : (v_final .- v_init) ./ v_scale

    # ---- u error ----
    p_u = plot(T, u_err,
               linewidth = 2,
               color      = :blue,
               xlabel     = "T",
               ylabel     = "Relative error",
               legend     = false,
               dpi        = dpi)
    save_plot(p_u, prefix, "_u_change.png", outdir)

    # ---- v error ----
    p_v = plot(T, v_err,
               linewidth = 2,
               color      = :red,
               xlabel     = "T",
               ylabel     = "Relative error",
               legend     = false,
               dpi        = dpi)
    save_plot(p_v, prefix, "_v_change.png", outdir)

    return p_u, p_v
end
evaluate_field(field::Function, T_grid, z_grid) = [field(zi, tj) for tj in T_grid, zi in z_grid]

function evaluate_field(field::AbstractMatrix, T_grid, z_grid)
    size(field) == (length(T_grid), length(z_grid)) ||
        throw(DimensionMismatch("field must have size $((length(T_grid), length(z_grid)))"))
    return field
end

function max_percentage_error_over_z(actual, reference, T_grid, z_grid)
    A = evaluate_field(actual, T_grid, z_grid)
    R = evaluate_field(reference, T_grid, z_grid)
    scale = maximum(abs, R)
    scale > 0 || throw(ArgumentError("reference is identically zero, percentage error is undefined"))
    return [100 * maximum(abs.(abs.(A[:, i]) .- abs.(R[:, i]))) / scale for i in eachindex(z_grid)]
end

"""
    plot_true_errors(; T_grid, z_grid, actual_u=nothing, reference_u=nothing,
                     actual_v=nothing, reference_v=nothing,
                     prefix="error plots", dpi=300, outdir="results")

    Plot the maximum percentage error over T versus z, for each component with a pair.

    At each z:
        error[i] = 100 * max over T of abs(abs(actual) - abs(reference)) / max over (T, z) of abs(reference)

    - `T_grid` : retarded time array (length N_T)
    - `z_grid` : propagation distances (length N_z)
    - `actual_u`, `reference_u` : u solution and its reference. Each is a complex matrix (N_T, N_z)
      or a callable `field(z, T)`. A component is plotted only when both of its entries are given
    - `actual_v`, `reference_v` : same for v

    If no pair is given, nothing is plotted. Returns the vector of plots.
"""
function plot_true_errors(;
    T_grid, z_grid,
    actual_u=nothing, reference_u=nothing,
    actual_v=nothing, reference_v=nothing,
    prefix="error plots", dpi=300, outdir="results"
)
    plots = []
    for (component, actual, reference, color) in (
        ("u", actual_u, reference_u, :blue),
        ("v", actual_v, reference_v, :red),
    )
        (isnothing(actual) || isnothing(reference)) && continue
        err = max_percentage_error_over_z(actual, reference, T_grid, z_grid)
        p = plot(
            z_grid, err;
            xlabel = "z",
            ylabel = "relative error, %",
            linewidth = 2,
            color  = color,
            legend = false,
            dpi    = dpi
        )
        save_plot(p, prefix, "_$(component)_error_vs_z.png", outdir)
        push!(plots, p)
    end
    return plots
end

"""
    fit_order(h, err)

Least-squares slope of `log10(err)` versus `log10(h)`. For an order-p scheme `err ∝ h^p`, so the slope is p.
"""
function fit_order(h, err)
    A = [ones(length(h)) log10.(h)]
    coef = A \ log10.(err)
    return coef[2]
end

function sci_label(x)
    e = floor(Int, log10(x))
    m = round(x / 10.0^e, digits=2)
    mantissa = m == round(m) ? string(Int(round(m))) : string(m)
    return "$(mantissa)e$(e)"
end

"""
    plot_convergence(h, err; prefix="validation_convergence", dpi=300, outdir="results")

    Log-log plot of the error against the step parameter `h`. The measurements are black
    triangles, and the least-squares line `log10(err) = k*log10(h) + c` is a black dotted line.
    The x ticks are the measured `h` values.

    - `h` : step parameters
    - `err` : errors for each `h`

    Returns the plot and the fitted slope `k`.
"""
function plot_convergence(h, err; xlabel="h", show_fit=true, prefix="validation_convergence", dpi=300, outdir="results")
    order = fit_order(h, err)
    intercept = sum(log10.(err) .- order .* log10.(h)) / length(h)
    p = plot(
        h, err;
        seriestype = :scatter,
        marker = :utriangle,
        markersize = 7,
        markercolor = :black,
        markerstrokecolor = :black,
        label = "measurements",
        xscale = :log10,
        yscale = :log10,
        xticks = (h, sci_label.(h)),
        xlabel = xlabel,
        ylabel = "max relative error, %",
        legend = :topleft,
        dpi    = dpi
    )
    if show_fit
        h_line = range(minimum(h), maximum(h), length=100)
        plot!(
            p, h_line, 10 .^(intercept .+ order .* log10.(h_line));
            label = "linear fit (k=$(round(order, digits=2)), c=$(round(intercept, digits=2)))",
            linestyle = :dot,
            linewidth = 2,
            color = :black,
        )
    end
    save_plot(p, prefix, ".png", outdir)
    return p, order
end

"""
    plot_energy_exchange(runs; prefix="energy_exchange", dpi=300, outdir="results")

    Plot the energies `E_u(z)` and `E_v(z)` of each run in a grid of panels, one panel per run.

    - `runs` : vector of named tuples `(label, z, E_u, E_v)`
"""
function plot_energy_exchange(runs; prefix="energy_exchange", dpi=300, outdir="results")
    panels = [
        begin
            p = plot(
                run.z, run.E_u;
                label = "E_u",
                xlabel = "z",
                ylabel = "energy",
                title = run.label,
                linewidth = 2,
                color = :blue,
                legend = :best,
            )
            plot!(p, run.z, run.E_v; label = "E_v", linewidth = 2, color = :red)
            p
        end
        for run in runs
    ]
    ncols = min(length(panels), 2)
    nrows = cld(length(panels), ncols)
    p = plot(panels...; layout = (nrows, ncols), size = (900, 400 * nrows), dpi = dpi)
    save_plot(p, prefix, ".png", outdir)
    return p
end

"""
    plot_energy_evolution(T, z, U, V; prefix="energy_evolution", dpi=300, outdir="results")

    Plot the energies `E_u(z) = ∫|u|² dt` and `E_v(z) = ∫|v|² dt` of the solution matrices
    `U`, `V` of size `(length(T), length(z))`.
"""
function plot_energy_evolution(T, z, U, V; prefix="energy_evolution", dpi=300, outdir="results")
    dT = T[2] - T[1]
    E_u = vec(dT .* sum(abs2, U; dims=1))
    E_v = vec(dT .* sum(abs2, V; dims=1))
    return plot_energy_exchange([(label = "u and v", z = z, E_u = E_u, E_v = E_v)]; prefix, dpi, outdir)
end

"""
    plot_energy_total(T, z, U, V; prefix="energy_total", dpi=300, outdir="results")

    Plot the total energy `E_u(z) + E_v(z)` of the solution matrices `U`, `V` of size
    `(length(T), length(z))`.
"""
function plot_energy_total(T, z, U, V; prefix="energy_total", dpi=300, outdir="results")
    dT = T[2] - T[1]
    E = vec(dT .* sum(abs2, U; dims=1)) .+ vec(dT .* sum(abs2, V; dims=1))
    p = plot(
        z, E;
        xlabel = "z",
        ylabel = "energy",
        title = "total energy, u and v",
        linewidth = 2,
        color = :black,
        legend = false,
        dpi = dpi,
    )
    save_plot(p, prefix, ".png", outdir)
    return p
end

fixed_ylims(::Nothing) = nothing
fixed_ylims(top::Real) = (0.0, float(top))
fixed_ylims(lims) = (float(lims[1]), float(lims[2]))

"""
    envelope_clip_recorder(Tspan, N_T; every=100, fps=20, dpi=100, ylims=nothing)

    Record the envelopes `|u(t)|` and `|v(t)|` during `solve_coupled_gl`. Returns a named tuple:
    - `on_step` : callback for the solver's `on_step` keyword. Every `every` steps it adds a frame.
    - `save`    : `save(path)` writes the recorded frames as an animation (e.g. `"envelope.gif"`).
    The time grid is the one of the solver. The y-axis is fixed for the whole clip: `ylims` if given,
    otherwise `(0, 1.2 × max)` of the first recorded frame.
"""
function envelope_clip_recorder(Tspan, N_T; every=100, fps=20, dpi=100, ylims=nothing)
    dT = (Tspan[2] - Tspan[1]) / N_T
    T = range(Tspan[1], step=dT, length=N_T)
    anim = Animation()
    y = Ref{Any}(fixed_ylims(ylims))
    function on_step(k, z_k, u, v)
        (k - 1) % every == 0 || return
        isnothing(y[]) && (y[] = (0.0, 1.2 * max(maximum(abs, u), maximum(abs, v))))
        p = plot(
            T, abs.(u);
            label = "|u|",
            xlabel = "t",
            ylabel = "envelope",
            title = "z = $(round(z_k, digits=2))",
            ylims = y[],
            linewidth = 2,
            color = :blue,
            dpi = dpi,
        )
        plot!(p, T, abs.(v); label = "|v|", linewidth = 2, color = :red)
        frame(anim, p)
    end
    save(path) = gif(anim, path; fps=fps, show_msg=false)
    return (; on_step, save)
end

"""
    plot_gamma_grid(runs; prefix="gamma_grid", cmap=:plasma, dpi=300, outdir="results", max_z_points=2000)

    Grid of panels with one row per run and three columns: the energies `E_u(z)`, `E_v(z)`
    (left), `|u(z,t)|` (middle) and `|v(z,t)|` (right).

    - `runs` : vector of named tuples `(label, T, z, E_u, E_v, U, V)`, with `U`, `V` of size
      `(length(T), length(z))`
"""
function plot_gamma_grid(runs; prefix="gamma_grid", cmap=:plasma, dpi=300, outdir="results", max_z_points=2000)
    panels = []
    for run in runs
        cols = 1:cld(length(run.z), max_z_points):length(run.z)
        z_plot = run.z[cols]
        pe = plot(
            run.z, run.E_u;
            label = "E_u", color = :blue, linewidth = 2,
            xlabel = "z", ylabel = "energy", title = run.label, legend = :best,
        )
        plot!(pe, run.z, run.E_v; label = "E_v", color = :red, linewidth = 2)
        pu = heatmap(z_plot, run.T, abs.(run.U[:, cols]); color = cmap, xlabel = "z", ylabel = "t", title = "|u|")
        pv = heatmap(z_plot, run.T, abs.(run.V[:, cols]); color = cmap, xlabel = "z", ylabel = "t", title = "|v|")
        append!(panels, [pe, pu, pv])
    end
    p = plot(panels...; layout = (length(runs), 3), size = (1500, 350 * length(runs)), dpi = dpi)
    save_plot(p, prefix, ".png", outdir)
    return p
end

end # module

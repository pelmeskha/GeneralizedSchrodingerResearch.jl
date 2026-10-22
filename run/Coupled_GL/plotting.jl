module CoupledGLPlotting

using Plots

export plot_coupled_gl, plot_profiles, plot_error_profiles, plot_true_errors

const RESULTS_DIR = normpath(joinpath(@__DIR__, "../..", "results", "Coupled_GL"))

function save_plot(p, prefix, suffix, outdir)
    mkpath(outdir)
    path = joinpath(outdir, string(prefix, suffix))
    savefig(p, path)
    return path
end

"""
    plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis, outdir=RESULTS_DIR)

    Save two heatmaps showing the evolution of |u| and |v| as functions of z and T.

    - `T` : vector of retarded time values (length N_T)
    - `z` : vector of propagation distances (length N_z)
    - `U` : complex array of size (N_T, N_z) containing u(z,T)
    - `V` : complex array of size (N_T, N_z) containing v(z,T)
    - `prefix` : string prepended to output filenames
    - `outdir` : directory for saved figures (created if missing)
    - `cmap`  : colormap (symbol), e.g., :plasma, :inferno, :turbo
"""
function plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis, outdir=RESULTS_DIR)
    # |u|
    p1 = heatmap(
        z,
        T,
        abs.(U),
        xlabel = "z", ylabel = "T",
        title  = "|u(z,T)|",
        color  = cmap,
        dpi    = 300
    )
    save_plot(p1, prefix, "_u.png", outdir)

    # |v|
    p2 = heatmap(
        z,
        T,
        abs.(V),
        xlabel = "z", ylabel = "T",
        title  = "|v(z,T)|",
        color  = cmap,
        dpi    = 300
    )
    save_plot(p2, prefix, "_v.png", outdir)
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
                  prefixes=nothing, prefix="profiles", dpi=300, outdir=RESULTS_DIR)

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
    - `prefix` : filename prefix for saved figures
    - `dpi` : resolution
    - `outdir` : directory for saved figures (created if missing)

    Returns the vector of plots.
"""
function plot_profiles(;
    T_grid, z_grid, z_timestamps,
    u_solutions=nothing, v_solutions=nothing,
    prefixes=nothing,
    prefix="profiles", dpi=300, outdir=RESULTS_DIR
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
                xlabel = "T",
                ylabel = "|$component|",
                title  = "$component at z=$z_label",
                legend = :best,
                dpi    = dpi
            )
            for (i, solution) in enumerate(solutions)
                profile = abs.(sample_profile(solution, T_grid, z_grid, z_value))
                plot!(
                    p, T_grid, profile;
                    label = "$component, $(names[i]), at z=$z_label",
                    linewidth = i == 1 ? 1.5 : 3.5,
                    linestyle = i == 1 ? :solid : :dot,
                    color = PROFILE_COLORS[mod1(i, length(PROFILE_COLORS))]
                )
            end
            save_plot(p, prefix, "_$(component)_z=$z_label.png", outdir)
            push!(plots, p)
        end
    end

    return plots
end
"""
    plot_error_profiles(T, z, U, V; prefix="error", dpi=300, outdir=RESULTS_DIR)

    Compute and plot the relative change of |u| and |v| from start to end.

    The relative error is defined as
        ( |field(T, z_end)| - |field(T, z_start)| ) / max_T |field(T, z_start)|
    so that the curve shows the change normalised to the initial peak amplitude.
"""
function plot_error_profiles(T, z, U, V; prefix="error", dpi=300, outdir=RESULTS_DIR)
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
               title      = "|u|: (final - initial) / max(|u₀|)",
               legend     = false,
               dpi        = dpi)
    save_plot(p_u, prefix, "_u_error.png", outdir)

    # ---- v error ----
    p_v = plot(T, v_err,
               linewidth = 2,
               color      = :red,
               xlabel     = "T",
               ylabel     = "Relative error",
               title      = "|v|: (final - initial) / max(|v₀|)",
               legend     = false,
               dpi        = dpi)
    save_plot(p_v, prefix, "_v_error.png", outdir)

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
                     prefix="error plots", dpi=300, outdir=RESULTS_DIR)

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
    prefix="error plots", dpi=300, outdir=RESULTS_DIR
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
            title  = "Max percentage error over z",
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

end # module

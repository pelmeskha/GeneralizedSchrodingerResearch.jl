using GeneralizedSchrodingerResearch.Solvers: solve_coupled_gl
using Plots

"""
    plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis)

    Save two heatmaps showing the evolution of |u| and |v| as functions of z and T.

    - `T` : vector of retarded time values (length N_T)
    - `z` : vector of propagation distances (length N_z)
    - `U` : complex array of size (N_T, N_z) containing u(z,T)
    - `V` : complex array of size (N_T, N_z) containing v(z,T)
    - `prefix` : string prepended to output filenames
    - `cmap`  : colormap (symbol), e.g., :plasma, :inferno, :turbo
"""
function plot_coupled_gl(T, z, U, V; prefix="coupled_gl", cmap=:viridis)
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
    savefig(p1, "$(prefix)_u.png")

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
    savefig(p2, "$(prefix)_v.png")
end
"""
    plot_profiles(T, z, U, V; prefix="profiles", dpi=300)

    Plot the amplitude profiles |u(T)| and |v(T)| at the start (z[1]) and end (z[end]).

    - `T` : retarded time array (length N_T)
    - `z` : propagation distances (length N_z)
    - `U` : complex array (N_T, N_z) for u
    - `V` : complex array (N_T, N_z) for v
    - `prefix` : filename prefix for saved figures
    - `dpi` : resolution
"""
function plot_profiles(T, z, U, V; prefix="profiles", dpi=300)
    # Initial and final indices
    i_start = 1
    i_end   = length(z)

    # Extract profiles
    u_init = abs.(U[:, i_start])
    u_final = abs.(U[:, i_end])
    v_init = abs.(V[:, i_start])
    v_final = abs.(V[:, i_end])

    # ---- u ----
    p_u = plot(T, u_init,
               label = "z = $(z[i_start])",
               linewidth = 2,
               xlabel = "T",
               ylabel = "|u|",
               title  = "u - initial vs final",
               legend = :best,
               dpi    = dpi)
    plot!(p_u, T, u_final,
          label = "z = $(z[i_end])",
          linewidth = 2)
    savefig(p_u, "$(prefix)_u_profiles.png")

    # ---- v ----
    p_v = plot(T, v_init,
               label = "z = $(z[i_start])",
               linewidth = 2,
               xlabel = "T",
               ylabel = "|v|",
               title  = "v - initial vs final",
               legend = :best,
               dpi    = dpi)
    plot!(p_v, T, v_final,
          label = "z = $(z[i_end])",
          linewidth = 2)
    savefig(p_v, "$(prefix)_v_profiles.png")

    return p_u, p_v
end
"""
    plot_error_profiles(T, z, U, V; prefix="error", dpi=300)

    Compute and plot the relative change of |u| and |v| from start to end.

    The relative error is defined as
        ( |field(T, z_end)| - |field(T, z_start)| ) / max_T |field(T, z_start)|
    so that the curve shows the change normalised to the initial peak amplitude.
"""
function plot_error_profiles(T, z, U, V; prefix="error", dpi=300)
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
    savefig(p_u, "$(prefix)_u_error.png")

    # ---- v error ----
    p_v = plot(T, v_err,
               linewidth = 2,
               color      = :red,
               xlabel     = "T",
               ylabel     = "Relative error",
               title      = "|v|: (final - initial) / max(|v₀|)",
               legend     = false,
               dpi        = dpi)
    savefig(p_v, "$(prefix)_v_error.png")

    return p_u, p_v
end

# pulse parameters
k = 0
c0 = 1
mu = 1
I1 = 0.1
ω = 0
u0(T) = 2*sqrt(2*I1)/abs(mu) * abs(1-2/(1+ exp(mu*(0-c0*T-0)))) * exp(1im * (ω*T))
v0(T) = 0 # 2*sqrt(2*I1)/abs(mu) - 2*sqrt(2*I1)/abs(mu) * abs(1-2/(1+ exp(mu*(0-c0*T-0)))) * exp(1im * (ω*T))

# equation parameters. Construct corresponding equation
β2 = -1.0; # Any
Ω_g = 20.0; # not matter

Δβ = (β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω)/(2*c0); # derieved from restrictions
γ = -mu^4*c0^2*β2/(32*I1) # derieved from restrictions

δ = (1 - β2*c0*ω / c0); # from restrictions
g = 0; # from restrictions


# domain
Tspan = (-10.0, 10.0)
zspan = (0.0, 5.0)
dz    = 0.005
N_T   = 512

T, z, u, v = solve_coupled_gl(
    zspan,
    Tspan,
    dz,
    N_T,
    u0,
    v0,
    Δβ,
    δ,
    β2,
    g,
    Ω_g,
    γ
)
plot_coupled_gl(T, z, u, v; prefix="my_run", cmap=:plasma)
plot_profiles(T, z, u, v; prefix="my_run")
plot_error_profiles(T, z, u, v; prefix="my_run")
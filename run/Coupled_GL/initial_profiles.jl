include("plotting.jl")
using .CoupledGLPlotting
using GeneralizedSchrodingerResearch.AnalyticalSolutions.CoupledGLAnalytical

LINE_COLOR = :black
LINE_WIDTH = 2.5
LINE_STYLE = :solid
SHOW_LEGEND = false

c0 = 1.4
k = 12.6
β2 = 0.027168436376788473
Ω_g = 66.08752874553595
mu = 0.97
I1 = 0.1

# tanh (dark soliton)
Tspan_tanh = (-10.0, 10.0)
tanh_args = ModelingArgs(
    (c0=c0, k=k, ω=0, β2=β2, Ω_g=Ω_g, w0=0, mu=mu, I1=I1),
    construct_dependent_args_for_tanh,
)

# sn (dark elliptic soliton), Δβ and γ from the same tanh restrictions
Tspan_sn = (-10.0, 10.0)
sn_shift = (Tspan_sn[2] - Tspan_sn[1]) / 4
sn_ω = 0
sn_mu = 1.3
sn_Δβ = (-β2*c0^3*sn_mu^2 + 2*β2*c0*sn_ω^2 + 4*c0*k - 4*sn_ω) / (2*c0)
sn_γ = sn_mu^4 * c0^2 * β2 / (32*I1)
sn_args = ModelingArgs(
    (c0=c0, k=k, ω=sn_ω, β2=β2, Ω_g=Ω_g, w0=0, Tspan=Tspan_sn, Δβ=sn_Δβ, γ=sn_γ),
    construct_dependent_args_for_sn,
)

# bright soliton (needs Δβ > 2k + β2ω² − 2ω/c0 and γ < 0)
Tspan_bright = (-10.0, 10.0)
bright_args = ModelingArgs(
    (c0=c0, k=k, ω=5, β2=β2, Ω_g=Ω_g, w0=0, Tspan=Tspan_bright, Δβ=19.327449097454128, γ=-0.013160409145296585),
    construct_dependent_args_for_bright,
)

cases = [
    ("tanh", Tspan_tanh, (z, t) -> tanh_solution(z, t, tanh_args)),
    ("sn", Tspan_sn, (z, t) -> sn_solution(z, t - sn_shift, sn_args)),
    ("bright", Tspan_bright, (z, t) -> bright_solution(z, t, bright_args)),
]

for (name, Tspan, u) in cases
    T_grid = collect(range(Tspan[1], Tspan[2], length=4001))
    plot_profiles(
        T_grid = T_grid,
        z_grid = [0.0],
        z_timestamps = [0.0],
        u_solutions = [u],
        prefixes = [name],
        colors = [LINE_COLOR],
        linewidths = [LINE_WIDTH],
        linestyles = [LINE_STYLE],
        legend = SHOW_LEGEND,
        prefix = "initial_$name",
    )
end

using GeneralizedSchrodingerResearch.Solvers: solve_coupled_gl
include("plotting.jl")
using .CoupledGLPlotting
using GeneralizedSchrodingerResearch.AnalyticalSolutions.CoupledGLAnalytical
using GeneralizedSchrodingerResearch.AnalyticalSolutions.CoupledGLResiduals
include("../../src/utilities/periodic_time.jl")

# domain
Tspan = (-10.0, 10.0)
zspan = (0.0, 50.0)
dz    = 0.0022231407344530316 # must be lower than dt. Scheme is 2 order on z
dt    = 20.0 / 950
N_T   = round(Int, (Tspan[2]-Tspan[1]) / dt)

# pulse parameters
k = 12.6
c0 = 1.4
mu = 0.97
I1 = 0.1
ω = 0 # sn needs ω·(Tspan length) ∈ 2πℤ for a periodic carrier; bright used ω = 5
w0 = 0 # displacement

# equation parameters. Construct corresponding equation
β2 = 0.027168436376788473; # should be ~ 0.025
Ω_g = 66.08752874553595; # should be ~ [10, 100]

# sn (dark) solution: Δβ and γ follow from the tanh restrictions ξ₁ = μ²/2, ξ₂ = −μ⁴/(16 I₁).
Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω) / (2*c0)
γ  = mu^4 * c0^2 * β2 / (32*I1)

args = ModelingArgs(
    (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
    construct_dependent_args_for_sn,
)

u0(T) = sn_solution(0.0, T, args)
v0(T) = 0

true_u(z, T) = sn_solution(z, T, args)
true_v(z, T) = 0

# Bright soliton (previous setup, with ω = 5):
# Δβ = 19.327449097454128
# γ  = -0.013160409145296585
# args = ModelingArgs(
#     (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
#     construct_dependent_args_for_bright,
# )
# u0(T) = bright_solution(0.0, T, args)
# true_u(z, T) = bright_solution(z, T, args)
# true_v(z, T) = 0

# Periodic analytical reference.
true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, c0, w0))
true_v_periodic(z, T) = true_v(z, periodic_time(z, T, Tspan, c0, w0))


println("β2 = $(args.β2), Ω_g = $(args.Ω_g), Δβ = $(args.Δβ), γ = $(args.γ), δ = $(args.δ)")


println("Примерный критерий сходимости: $(dz/( (Tspan[2]-Tspan[1]) /N_T)^2) << $(0.2/(abs(β2)*pi^2)) ?")

T, z, u, v = solve_coupled_gl(zspan, Tspan, dz, N_T, u0, v0, args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ)

u_metrics = field_error(u, T, z, true_u_periodic)
println(
    "u vs analytical: max modulus error = ", round(u_metrics.max_modulus_error, sigdigits=3),
    ", max complex error = ", round(u_metrics.max_complex_error, sigdigits=3),
    ", max relative error = ", round(u_metrics.max_relative_error_percent, sigdigits=3), " %"
)

# Plotting is slow; uncomment only to inspect a run by eye.
plot_coupled_gl(T, z, u, v; prefix="2D plot actual", cmap=:plasma)
plot_coupled_gl(T, z, true_u_periodic, true_v_periodic; prefix="2D plot analytical", cmap=:plasma)
plot_profiles(
    T_grid = T,
    z_grid = z,
    z_timestamps = [zspan[1], zspan[2]],
    u_solutions = [u, true_u_periodic],
    # v_solutions = [v, true_v_periodic], # commented = empty array by default. No plots for v component
    prefixes=["actual", "analytical"],
)
plot_true_errors(
    T_grid = T,
    z_grid = z,
    actual_u = u,
    # actual_v = v,
    reference_u = true_u_periodic,
    # reference_v = true_v_periodic;
    prefix="error"
)
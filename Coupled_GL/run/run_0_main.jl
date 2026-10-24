using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, RESULTS_DIR
include(joinpath(@__DIR__, "case_setup.jl"))
using GeneralizedSchrodingerResearch.Utilities: field_error
using GeneralizedSchrodingerResearch.Plotting: plot_coupled_gl, plot_profiles, plot_true_errors

setup = reference()
(; Tspan, zspan, args, u0, v0, true_u_periodic, true_v_periodic) = setup

dz    = 0.0022 # must be lower than dt. Scheme is 2 order on z
dt    = 0.021
N_T   = round(Int, (Tspan[2]-Tspan[1]) / dt)

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

println("β2 = $(args.β2), Ω_g = $(args.Ω_g), Δβ = $(args.Δβ), γ = $(args.γ), δ = $(args.δ)")
println("Примерный критерий сходимости: $(dz/( (Tspan[2]-Tspan[1]) /N_T)^2) << $(0.2/(abs(args.β2)*pi^2)) ?")

T, z, u, v = solve_coupled_gl(zspan, Tspan, dz, N_T, u0, v0, args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ)

true_u_filed = [true_u_periodic(zi, tj) for tj in T, zi in z]
true_v_filed = [true_v_periodic(zi, tj) for tj in T, zi in z]

u_metrics = field_error(u, T, z, true_u_filed)
println(
    "u vs analytical: max modulus error = ", round(u_metrics.max_modulus_error, sigdigits=3),
    ", max complex error = ", round(u_metrics.max_complex_error, sigdigits=3),
    ", max relative error = ", round(u_metrics.max_relative_error_percent, sigdigits=3), " %"
)

# Plotting is slow; uncomment only to inspect a run by eye.
plot_coupled_gl(T, z, u, v; prefix="2D plot actual", cmap=:plasma, outdir=RESULTS_DIR)
plot_coupled_gl(T, z, true_u_filed, true_v_filed; prefix="2D plot analytical", cmap=:plasma, outdir=RESULTS_DIR)
plot_profiles(
    T_grid = T,
    z_grid = z,
    z_timestamps = [zspan[1], zspan[2]],
    u_solutions = [u, true_u_filed],
    # v_solutions = [v, true_v_periodic], # commented = empty array by default. No plots for v component
    prefixes=["actual", "analytical"],
    outdir=RESULTS_DIR,
)
plot_true_errors(
    T_grid = T,
    z_grid = z,
    actual_u = u,
    # actual_v = v,
    reference_u = true_u_filed,
    # reference_v = true_v_periodic;
    prefix="error",
    outdir=RESULTS_DIR,
)

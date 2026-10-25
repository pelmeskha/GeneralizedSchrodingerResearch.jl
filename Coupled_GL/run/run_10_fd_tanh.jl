using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl_fd, RESULTS_DIR
include(joinpath(@__DIR__, "case_setup.jl"))
using GeneralizedSchrodingerResearch.Utilities: field_error
using GeneralizedSchrodingerResearch.Plotting: plot_coupled_gl, plot_profiles, plot_true_errors, plot_energy_evolution

setup = tanh_open_reference()
(; Tspan, zspan, args, u0, v0, true_u, true_v) = setup

BC    = :neumann   # the kink tends to ±A at the ends, so zero-flux (2nd kind) is the appropriate condition
dz    = 0.001      # finite-difference scheme is 2nd order in dz and in dt
dt    = 0.0125
N_T   = round(Int, (Tspan[2]-Tspan[1]) / dt) + 1
store_every = 10   # keep every 0.01 in z: the full fields would not fit in memory

println("bc = $BC, β2 = $(args.β2), Ω_g = $(args.Ω_g), Δβ = $(args.Δβ), γ = $(args.γ), δ = $(args.δ), μ = $(args.mu)")

T, z, u, v = solve_coupled_gl_fd(zspan, Tspan, dz, N_T, u0, v0, args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ; bc=BC, store_every=store_every)

true_u_filed = [true_u(zi, tj) for tj in T, zi in z]
true_v_filed = [true_v(zi, tj) for tj in T, zi in z]

u_metrics = field_error(u, T, z, true_u_filed)
println(
    "u vs analytical: max modulus error = ", round(u_metrics.max_modulus_error, sigdigits=3),
    ", max complex error = ", round(u_metrics.max_complex_error, sigdigits=3),
    ", max relative error = ", round(u_metrics.max_relative_error_percent, sigdigits=3), " %"
)

v_metrics = field_error(v, T, z, true_v_filed)
println(
    "v vs analytical: max modulus error = ", round(v_metrics.max_modulus_error, sigdigits=3),
    ", max complex error = ", round(v_metrics.max_complex_error, sigdigits=3),
    ", max relative error = ", round(v_metrics.max_relative_error_percent, sigdigits=3), " %"
)

# Plotting is slow; uncomment only to inspect a run by eye.
plot_coupled_gl(T, z, u, v; prefix="2D plot FD tanh $(BC)", cmap=:plasma, outdir=RESULTS_DIR)
plot_energy_evolution(T, z, u, v; prefix="energy FD tanh $(BC)", outdir=RESULTS_DIR)
plot_profiles(
    T_grid = T,
    z_grid = z,
    z_timestamps = [zspan[1], zspan[2]],
    u_solutions = [u, true_u_filed],
    v_solutions = [v, true_v_filed],
    prefixes=["FD tanh $(BC)", "analytical"],
    outdir=RESULTS_DIR,
)
plot_true_errors(
    T_grid = T,
    z_grid = z,
    actual_u = u,
    reference_u = true_u_filed;
    prefix="error FD tanh $(BC)",
    outdir=RESULTS_DIR,
)

using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, RESULTS_DIR
include(joinpath(@__DIR__, "case_setup.jl"))
using GeneralizedSchrodingerResearch.Utilities: field_error
using GeneralizedSchrodingerResearch.Plotting: plot_coupled_gl, plot_profiles, plot_true_errors, plot_energy_evolution, plot_energy_total, envelope_clip_recorder

L       = 5.0
Tspan = (-L, L)
zspan = (0.0, 1000.0)
dz    = 0.0022 # must be lower than dt. Scheme is 2 order on z
dt    = 0.021
N_T   = round(Int, (Tspan[2]-Tspan[1]) / dt)

# Pulse parameters
c0    = 1.5
k     = 12.8
ω     = 4.2
β2    = 0.04
Ω_g   = 10
w0    = 0
Δβ    = 21
γ     = -0.01

args = ModelingArgs(
    (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
    construct_dependent_args_for_bright,
)

u0(T) = bright_solution(0.0, T, args)
v0(T) = 0.1 * exp(1im * pi/4) * bright_solution(0.0, T, args)

# Equation parameters
Δβ = args.Δβ
δ = args.δ - 1/args.c0
β2 = args.β2
g = 0.0002 # args.g
Ω_g = args.Ω_g
γ = -0.05 # [0, -0.01, -0.02, -0.04]

clip = envelope_clip_recorder(Tspan, N_T; every=3000, ylims=25)
T, z, u, v = solve_coupled_gl(
    zspan, Tspan, dz, N_T, u0, v0, Δβ, δ, β2, g, Ω_g, γ;
    on_step=clip.on_step
)
clip.save(joinpath(RESULTS_DIR, "envelope.gif"))

plot_energy_evolution(T, z, u, v; prefix="energy", outdir=RESULTS_DIR)
# plot_energy_total(T, z, u, v; prefix="energy_total", outdir=RESULTS_DIR)

# Plotting is slow; uncomment only to inspect a run by eye.
plot_coupled_gl(T, z, u, v; prefix="2D plot actual", cmap=:plasma, outdir=RESULTS_DIR)
plot_profiles(
    T_grid = T,
    z_grid = z,
    z_timestamps = [zspan[1], zspan[2]],
    u_solutions = [u],
    v_solutions = [v],
    prefixes=["actual"],
    outdir=RESULTS_DIR,
)

using GeneralizedSchrodingerResearch.CoupledGL: max_relative_error_over_run_fd, save_convergence_csv, RESULTS_DIR
include(joinpath(@__DIR__, "case_setup.jl"))
using GeneralizedSchrodingerResearch.Plotting: plot_convergence
using Printf

# Convergence of the finite-difference solver. Expected: 2nd order in dz and in dT for both boundary kinds.
#
# Study A (exact linear eigenmodes, both boundary kinds): with δ = 0 the operator has no first-derivative
# term, so u = exp(λz) cos(kπ(T-Tmin)/L) (Neumann) and u = exp(λz) sin(kπ(T-Tmin)/L) (Dirichlet) are exact.
# The Neumann mode is nonzero at the ends, so it tests the boundary treatment directly.
#
# Study B (soliton): the circular vector soliton on the open domain, dz fixed small.

const MODE_ARGS = (Δβ=0.7, δ=0.0, β2=0.3, g=0.4, Ω_g=1.5, γ=0.0)
const MODE_TSPAN = (-5.0, 5.0)
const MODE_ZSPAN = (0.0, 1.0)
const MODE_K = 1

function mode_setup(bc)
    Tmin, Tmax = MODE_TSPAN
    L = Tmax - Tmin
    kk = MODE_K * π / L
    a2 = MODE_ARGS.g/(2*MODE_ARGS.Ω_g^2) - 1im*MODE_ARGS.β2/2
    λ  = MODE_ARGS.g/2 + 1im*MODE_ARGS.Δβ/2 - a2*kk^2
    spatial(T) = bc === :neumann ? cos(kk*(T - Tmin)) : sin(kk*(T - Tmin))
    true_u(z, T) = exp(λ*z) * spatial(T)
    true_v(z, T) = 0
    return (; tag="mode_$(bc)", Tspan=MODE_TSPAN, zspan=MODE_ZSPAN, args=MODE_ARGS,
        u0=T -> spatial(T), v0=T -> 0.0, true_u, true_v)
end

function soliton_setup()
    setup = vector_circular_open_reference()
    return setup
end

function run_study(setup, bc; vary, values, fixed, xlabel, tag, expected_order)
    rows = []
    for x in values
        dz = vary === :dz ? x : fixed.dz
        dt = vary === :dt ? x : fixed.dt
        N_T = round(Int, (setup.Tspan[2] - setup.Tspan[1]) / dt) + 1
        r = max_relative_error_over_run_fd(setup, dz, N_T; bc=bc)
        @printf("bc = %s, dz = %.5f, dT = %.5f, N_T = %d, N_z = %d, max relative error = %.6f %% (%.1f s)\n",
                bc, dz, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime)
        push!(rows, [dz, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime])
    end
    save_convergence_csv(
        joinpath(RESULTS_DIR, "fd_convergence_$(tag)_$(bc).csv"),
        "dz,dT,N_T,N_z,max_relative_error_percent,runtime_s",
        rows,
    )
    x_values = vary === :dz ? [r[1] for r in rows] : [r[2] for r in rows]
    _, order = plot_convergence(
        x_values, [r[5] for r in rows];
        xlabel=xlabel, prefix="fd_convergence_$(tag)_$(bc)", outdir=RESULTS_DIR,
    )
    println("$tag study, bc = $bc: fitted order = $(round(order, digits=3)) (expected $expected_order)")
    return order
end

for bc in (:dirichlet, :neumann)
    setup = mode_setup(bc)
    run_study(setup, bc; vary=:dz, values=[0.08, 0.04, 0.02, 0.01], fixed=(dt=0.0125, dz=0.0),
        xlabel="dz", tag="mode_dz", expected_order=2)
    run_study(setup, bc; vary=:dt, values=[0.4, 0.2, 0.1, 0.05, 0.025], fixed=(dt=0.0, dz=0.001),
        xlabel="dT", tag="mode_dt", expected_order=2)
end

soliton = soliton_setup()
run_study(soliton, :dirichlet; vary=:dt, values=[0.05, 0.025, 0.0125], fixed=(dt=0.0, dz=0.0005),
    xlabel="dT", tag="soliton_dt", expected_order=2)
println("Results saved to $RESULTS_DIR")

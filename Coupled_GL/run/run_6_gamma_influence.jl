using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, RESULTS_DIR
using GeneralizedSchrodingerResearch.Plotting: plot_energy_exchange
include(joinpath(@__DIR__, "case_setup.jl"))

setup = bright_reference()
(; Tspan, zspan, args, u0, v0) = setup

dz = 0.0022231407344530316
dt = 20.0 / 950
N_T = round(Int, (Tspan[2]-Tspan[1]) / dt)
dT = (Tspan[2]-Tspan[1]) / N_T

γ_reference = args.γ
γ_values = [0.0, γ_reference / 2, γ_reference, 2 * γ_reference]

seed = 0.1
seed_phase = π / 4
v0_seeded(T) = seed * cis(seed_phase) * u0(T)

function energy_history(γ)
    E_u = Float64[]
    E_v = Float64[]
    record!(k, z_k, u, v) = (push!(E_u, dT * sum(abs2, u)); push!(E_v, dT * sum(abs2, v)))
    _, z, _, _ = solve_coupled_gl(
        zspan, Tspan, dz, N_T, u0, v0_seeded,
        args.Δβ, args.δ, args.β2, args.g, args.Ω_g, γ;
        show_progress=false, store=false, on_step=record!,
    )
    return z, E_u, E_v
end

runs = map(γ_values) do γ
    z, E_u, E_v = energy_history(γ)
    println("γ = $γ: E_u(0) = $(round(E_u[1], sigdigits=4)), E_u(end) = $(round(E_u[end], sigdigits=4)), max E_v = $(round(maximum(E_v), sigdigits=4))")
    (label = "γ = $(round(γ, sigdigits=3))", z = z, E_u = E_u, E_v = E_v)
end

plot_energy_exchange(runs; prefix="gamma_influence", outdir=RESULTS_DIR)

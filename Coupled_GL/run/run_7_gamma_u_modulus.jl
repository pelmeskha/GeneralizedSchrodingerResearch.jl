using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, RESULTS_DIR
using GeneralizedSchrodingerResearch.Plotting: plot_coupled_gl
include(joinpath(@__DIR__, "case_setup.jl"))

L     = 6.0
Tspan = (-L, L)
zspan = (0.0, 500.0)
dz    = 0.0022
dt    = 0.021
N_T   = round(Int, (Tspan[2]-Tspan[1]) / dt)
dT    = (Tspan[2]-Tspan[1]) / N_T
T     = collect(range(Tspan[1], step=dT, length=N_T))

c0    = 1.5
k     = 12.8
ω     = 4.2
β2    = 0.04
Ω_g   = 10
w0    = 0
Δβ    = 21
γ_ref = -0.01

args = ModelingArgs(
    (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ_ref),
    construct_dependent_args_for_bright,
)

u0(T) = bright_solution(0.0, T, args)
v0(T) = 0.1 * exp(1im * pi/4) * bright_solution(0.0, T, args)

δ = args.δ - 1/args.c0
g = args.g
γ_values = [-0.02, -0.04]
stride = 100

function run_gamma(γ)
    z_s = Float64[]
    U_cols = Vector{Vector{ComplexF64}}()
    V_cols = Vector{Vector{ComplexF64}}()
    function record!(k, z_k, u, v)
        (k - 1) % stride == 0 || return
        push!(z_s, z_k)
        push!(U_cols, copy(u))
        push!(V_cols, copy(v))
    end
    solve_coupled_gl(
        zspan, Tspan, dz, N_T, u0, v0, Δβ, δ, args.β2, g, Ω_g, γ;
        show_progress=false, store=false, on_step=record!,
    )
    return (z = z_s, U = reduce(hcat, U_cols), V = reduce(hcat, V_cols))
end

for γ in γ_values
    r = run_gamma(γ)
    plot_coupled_gl(T, r.z, r.U, r.V; prefix="2D plot gamma=$(γ) T=6", cmap=:plasma, outdir=RESULTS_DIR)
end

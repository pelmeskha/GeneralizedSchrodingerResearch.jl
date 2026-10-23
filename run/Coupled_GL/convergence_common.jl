using GeneralizedSchrodingerResearch.Solvers: solve_coupled_gl
include("plotting.jl")
using .CoupledGLPlotting
using GeneralizedSchrodingerResearch.AnalyticalSolutions.CoupledGLAnalytical
include("../../src/utilities/periodic_time.jl")

using Printf

const SOLUTION_TAG = "sn"

# domain
Tspan = (-10.0, 10.0)
zspan = (0.0, 50.0)

# pulse parameters
k = 12.6
c0 = 1.4
mu = 0.97
I1 = 0.1
ω = 0
w0 = 0

# equation parameters
β2 = 0.027168436376788473
Ω_g = 66.08752874553595
Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω) / (2*c0)
γ  = mu^4 * c0^2 * β2 / (32*I1)

args = ModelingArgs(
    (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
    construct_dependent_args_for_sn,
)

u0(T) = sn_solution(0.0, T, args)
v0(T) = 0

true_u(z, T) = sn_solution(z, T, args)
true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, c0, w0))

const RESULTS_CONV_DIR = normpath(joinpath(@__DIR__, "../../results/Coupled_GL"))

function max_relative_error_over_run(zspan, Tspan, dz, N_T)
    dT = (Tspan[2] - Tspan[1]) / N_T
    T = collect(range(Tspan[1], step=dT, length=N_T))
    max_err = 0.0
    max_ref = 0.0
    check_every = max(1, round(Int, 0.01 / dz))

    function check_step!(k, z_k, u, v)
        (k - 1) % check_every == 0 || return
        for j in eachindex(T)
            ref = true_u_periodic(z_k, T[j])
            max_err = max(max_err, abs(u[j] - ref))
            max_ref = max(max_ref, abs(ref))
        end
    end

    t0 = time()
    _, z, _, _ = solve_coupled_gl(
        zspan, Tspan, dz, N_T, u0, v0,
        args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ;
        show_progress=false, store=false, on_step=check_step!,
    )
    runtime = time() - t0

    return (
        dT=dT,
        N_z=length(z),
        max_relative_error_percent=100 * max_err / max_ref,
        runtime=runtime,
    )
end

function save_convergence_csv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for r in rows
            println(io, join(r, ","))
        end
    end
end

using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, RESULTS_DIR
using GeneralizedSchrodingerResearch.Utilities: field_error, periodic_time
using GeneralizedSchrodingerResearch.CoupledGL.CoupledGLAnalytical
using Random

const N_ROUNDS = 20
const SEED = 20261022

const Tspan = (-10.0, 10.0)
const zspan = (0.0, 50.0)

const dz = 0.0022231407344530316
const N_T = 950

const k = 12.6
const c0 = 1.4
const ω = 5
const w0 = 0

const BETA2_RANGE = (0.01, 0.05)
const OMEGA_G_RANGE = (10.0, 100.0)
const GAMMA_ABS_RANGE = (0.001, 0.05)
const DELTA_BETA_MAX = 40.0

const CSV_PATH = joinpath(RESULTS_DIR, "random_search.csv")

log_uniform(rng, lo, hi) = exp(log(lo) + rand(rng) * (log(hi) - log(lo)))

function sample_args(rng)
    β2 = BETA2_RANGE[1] + rand(rng) * (BETA2_RANGE[2] - BETA2_RANGE[1])
    Ω_g = log_uniform(rng, OMEGA_G_RANGE...)
    γ = -log_uniform(rng, GAMMA_ABS_RANGE...)
    Δβ_min = 2k + β2 * ω^2 - 2ω / c0
    Δβ = Δβ_min + rand(rng) * (DELTA_BETA_MAX - Δβ_min)
    return ModelingArgs(
        (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
        construct_dependent_args_for_bright,
    )
end

function random_search()
    rng = MersenneTwister(SEED)

    mkpath(dirname(CSV_PATH))
    open(CSV_PATH, "w") do io
        println(io, "round,beta2,Omega_g,Delta_beta,gamma,delta,g,dz,N_T,dT,N_z,max_relative_error_percent,max_complex_error,max_modulus_error,runtime_s")
    end

    best_r, best_err = 0, Inf

    for r in 1:N_ROUNDS
        args = sample_args(rng)
        args.g == 0 || error("g must be 0, got $(args.g)")

        u0(T) = bright_solution(0.0, T, args)
        v0(T) = 0
        true_u_periodic(z, T) = bright_solution(z, periodic_time(z, T, Tspan, c0, w0), args)

        t0 = time()
        Tg, zg, uu, vv = solve_coupled_gl(
            zspan, Tspan, dz, N_T, u0, v0,
            args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ;
            show_progress=false,
        )
        m = field_error(uu, Tg, zg, true_u_periodic)
        runtime = time() - t0
        err = isnan(m.max_relative_error_percent) ? Inf : m.max_relative_error_percent

        row = [r, args.β2, args.Ω_g, args.Δβ, args.γ, args.δ, args.g,
               dz, N_T, Tg[2] - Tg[1], length(zg),
               m.max_relative_error_percent, m.max_complex_error, m.max_modulus_error, runtime]
        open(CSV_PATH, "a") do io
            println(io, join(row, ","))
        end

        if err < best_err
            best_r, best_err = r, err
        end
        println("round $r: β2 = $(round(args.β2, sigdigits=4)), Ω_g = $(round(args.Ω_g, sigdigits=4)), Δβ = $(round(args.Δβ, sigdigits=4)), γ = $(round(args.γ, sigdigits=4)), max rel. error = $(round(m.max_relative_error_percent, sigdigits=4)) % ($(round(runtime, digits=1)) s); best so far $(round(best_err, sigdigits=4)) % at round $best_r")

        uu = vv = nothing
        GC.gc()
    end

    println("Best: round $best_r, max relative error = $best_err %")
    println("History saved to $CSV_PATH")
end

random_search()

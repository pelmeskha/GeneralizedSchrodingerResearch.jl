using GeneralizedSchrodingerResearch.CoupledGL: solve_coupled_gl, sn_reference, RESULTS_DIR
using GeneralizedSchrodingerResearch.Utilities: field_error
using Random

const N_ROUNDS = 50
const DZ_RANGE = (0.002, 0.02)
const DT_RANGE = (0.02, 0.075)
const CSV_PATH = joinpath(RESULTS_DIR, "grid_search.csv")

setup = sn_reference()
(; Tspan, zspan, args, u0, v0, true_u_periodic) = setup

function grid_search()
    rng = MersenneTwister(20261021)
    sample_dz() = exp(log(DZ_RANGE[1]) + rand(rng) * (log(DZ_RANGE[2]) - log(DZ_RANGE[1])))
    sample_dt() = exp(log(DT_RANGE[1]) + rand(rng) * (log(DT_RANGE[2]) - log(DT_RANGE[1])))

    mkpath(dirname(CSV_PATH))
    open(CSV_PATH, "w") do io
        println(io, "round,dz,N_T,dT,N_z,max_relative_error_percent,max_complex_error,max_modulus_error,runtime_s")
    end

    best_r, best_dz, best_NT, best_err = 0, NaN, 0, Inf

    for r in 1:N_ROUNDS
        dz = sample_dz()
        dt = sample_dt()
        N_T = round(Int, (Tspan[2] - Tspan[1]) / dt)
        t0 = time()
        Tg, zg, uu, vv = solve_coupled_gl(
            zspan, Tspan, dz, N_T, u0, v0,
            args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ;
            show_progress=false,
        )
        m = field_error(uu, Tg, zg, true_u_periodic)
        runtime = time() - t0
        err = isnan(m.max_relative_error_percent) ? Inf : m.max_relative_error_percent

        row = [r, dz, N_T, Tg[2] - Tg[1], length(zg),
               m.max_relative_error_percent, m.max_complex_error, m.max_modulus_error, runtime]
        open(CSV_PATH, "a") do io
            println(io, join(row, ","))
        end

        if err < best_err
            best_r, best_dz, best_NT, best_err = r, dz, N_T, err
        end
        println("round $r: dz = $(round(dz, sigdigits=4)), N_T = $N_T, max rel. error = $(round(m.max_relative_error_percent, sigdigits=4)) % ($(round(runtime, digits=1)) s); best so far $(round(best_err, sigdigits=4)) % at dz = $(round(best_dz, sigdigits=4)), N_T = $best_NT")

        uu = vv = nothing
        GC.gc()
    end

    println("Best: round $best_r, dz = $best_dz, N_T = $best_NT, max relative error = $best_err %")
    println("History saved to $CSV_PATH")
end

grid_search()

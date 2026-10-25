using GeneralizedSchrodingerResearch.CoupledGL: max_relative_error_over_run_fd, save_convergence_csv, RESULTS_DIR
include(joinpath(@__DIR__, "case_setup.jl"))
using GeneralizedSchrodingerResearch.Plotting: plot_convergence
using Printf

# Grid convergence of the finite-difference solver for the tanh (dark soliton) solution with Neumann boundaries.
# Δz is fixed and small, Δt is varied. Expected: 2nd order in Δt.

const BC = :neumann
const DZ_FIXED = 0.00025
const DT_LIST = [0.05, 0.025, 0.0125]
const EXPECTED_ORDER = 2

setup = tanh_open_reference()
CSV_PATH = joinpath(RESULTS_DIR, "convergence_tanh_fd_$(BC).csv")

function run_convergence(setup, dt_list)
    rows = []
    for dt in dt_list
        N_T = round(Int, (setup.Tspan[2] - setup.Tspan[1]) / dt) + 1
        r = max_relative_error_over_run_fd(setup, DZ_FIXED, N_T; bc=BC)
        @printf("dt = %.4f, dz = %.5f, N_T = %d, N_z = %d, max relative error = %.6f %% (%.1f s)\n",
                r.dT, DZ_FIXED, N_T, r.N_z, r.max_relative_error_percent, r.runtime)
        push!(rows, [r.dT, DZ_FIXED, N_T, r.N_z, r.max_relative_error_percent, r.runtime])
    end
    return rows
end

rows = run_convergence(setup, DT_LIST)
save_convergence_csv(
    CSV_PATH,
    "dT,dz,N_T,N_z,max_relative_error_percent,runtime_s",
    rows,
)

dt_values = [r[1] for r in rows]
errors = [r[5] for r in rows]
_, order = plot_convergence(dt_values, errors; xlabel="dt", prefix="convergence_tanh_fd_$(BC)", outdir=RESULTS_DIR)
println("fitted order in dt = $(round(order, digits=3)) (expected $EXPECTED_ORDER)")
println("History saved to $CSV_PATH")

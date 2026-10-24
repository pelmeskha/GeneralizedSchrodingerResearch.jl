using GeneralizedSchrodingerResearch.CoupledGL: sn_reference, max_relative_error_over_run, save_convergence_csv, RESULTS_DIR
using GeneralizedSchrodingerResearch.Plotting: plot_convergence
using Printf

# Convergence study: the step parameter h sets both steps.
# Guess: dz = h^2 (Strang splitting is 2nd order in dz, so the error ~ dz^2 = h^4)
# and dt = h (the Fourier method in T is spectral, so dt only has to resolve the pulse).
const H_LIST = [0.2, 0.1, 0.05]
# const H_LIST = [0.02, 0.01, 0.005] # bright soliton; for sn the finer steps reach the round-off floor
dz_of(h) = h^2
dt_of(h) = h
const EXPECTED_ORDER = 4

setup = sn_reference()
CSV_PATH = joinpath(RESULTS_DIR, "validation_convergence_$(setup.tag).csv")

function run_convergence(setup, h_list)
    rows = []
    for h in h_list
        dz = dz_of(h)
        N_T = round(Int, (setup.Tspan[2] - setup.Tspan[1]) / dt_of(h))
        r = max_relative_error_over_run(setup, dz, N_T)
        push!(rows, (h=h, dz=dz, dT=r.dT, N_T=N_T, N_z=r.N_z,
                     err=r.max_relative_error_percent, runtime=r.runtime))
        @printf("h = %.4f, dz = %.3e, dt = %.4f, N_T = %d, N_z = %d, max relative error = %.5f %% (%.1f s)\n",
                h, dz, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime)
    end
    return rows
end

rows = run_convergence(setup, H_LIST)
save_convergence_csv(
    CSV_PATH,
    "h,dz,dT,N_T,N_z,max_relative_error_percent,runtime_s",
    [[r.h, r.dz, r.dT, r.N_T, r.N_z, r.err, r.runtime] for r in rows],
)

h_values = [r.h for r in rows]
errors = [r.err for r in rows]
p, order = plot_convergence(h_values, errors; prefix="validation_convergence_$(setup.tag)", outdir=RESULTS_DIR)
println("fitted order = $(round(order, digits=3)) (expected $EXPECTED_ORDER)")
println("History saved to $CSV_PATH")

include("convergence_common.jl")

# Study 1: dt fixed and well resolved, dz varied. Expected: error ~ dz^2.
const DT_FIXED = 0.05
const DZ_LIST = [0.04, 0.02, 0.01, 0.005, 0.0025]
# bright soliton: const DT_FIXED = 0.005; const DZ_LIST = [0.02, 0.01, 0.005, 0.0025, 0.00125]
const EXPECTED_ORDER_DZ = 2

# Study 2: dz fixed small, dt varied. Expected: error falls quickly and then flattens (spectral in t).
const DZ_FIXED = 0.0025
const DT_LIST = [0.2, 0.1, 0.05, 0.025]
# bright soliton: const DZ_FIXED = 1e-4; const DT_LIST = [0.04, 0.02, 0.01, 0.005]

function study_dz()
    N_T = round(Int, (Tspan[2] - Tspan[1]) / DT_FIXED)
    rows = []
    for dz in DZ_LIST
        r = max_relative_error_over_run(zspan, Tspan, dz, N_T)
        @printf("dz = %.5f, dt = %.4f, N_T = %d, N_z = %d, max relative error = %.6f %% (%.1f s)\n",
                dz, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime)
        push!(rows, [dz, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime])
    end
    save_convergence_csv(
        joinpath(RESULTS_CONV_DIR, "convergence_dz_$(SOLUTION_TAG).csv"),
        "dz,dT,N_T,N_z,max_relative_error_percent,runtime_s",
        rows,
    )
    _, order = plot_convergence(
        DZ_LIST, [r[5] for r in rows];
        xlabel="dz", prefix="convergence_dz_$(SOLUTION_TAG)", outdir=RESULTS_CONV_DIR,
    )
    println("dz study: fitted order in dz = $(round(order, digits=3)) (expected $EXPECTED_ORDER_DZ)")
end

function study_dt()
    rows = []
    for dt in DT_LIST
        N_T = round(Int, (Tspan[2] - Tspan[1]) / dt)
        r = max_relative_error_over_run(zspan, Tspan, DZ_FIXED, N_T)
        @printf("dt = %.4f, N_T = %d, N_z = %d, max relative error = %.6f %% (%.1f s)\n",
                r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime)
        push!(rows, [dt, r.dT, N_T, r.N_z, r.max_relative_error_percent, r.runtime])
    end
    save_convergence_csv(
        joinpath(RESULTS_CONV_DIR, "convergence_dt_$(SOLUTION_TAG).csv"),
        "dt,dT,N_T,N_z,max_relative_error_percent,runtime_s",
        rows,
    )
    plot_convergence(
        DT_LIST, [r[5] for r in rows];
        xlabel="dt", show_fit=false, prefix="convergence_dt_$(SOLUTION_TAG)", outdir=RESULTS_CONV_DIR,
    )
end

study_dz()
study_dt()
println("Results saved to $RESULTS_CONV_DIR")

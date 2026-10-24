module GNSE_3_5_7

using LinearAlgebra
using ProgressMeter
using Optim, Interpolations
using ..BaseUtilities: create_FD_matrix, integral_1, integral_2, filtration,
    find_threshold, shift_pulse_to_center, relative_error_to_amplitude, M0_M1_to_ε2_ε3, ε2_ε3_to_M0_M1

include("analytical_solutions/NSE/NSE_3_soliton.jl")
export NSE_3_soliton

include("analytical_solutions/NSE/NSE_3_5_soliton.jl")
export NSE_3_5_soliton

include("analytical_solutions/NSE/NSE_3_5_7_soliton.jl")
export precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton

include("solvers/solve.jl")
export solve

# Implemented in the CUDA package extension (ext/); available after `using CUDA`.
function cuda_solve end
export cuda_solve

# Implemented in the GLMakie package extension (ext/GLMakie/); available after `using GLMakie`.
function live_plot_record end
export live_plot_record

include("utilities/construct_approximate_solutions.jl")
export construct_approximate_NSE_3_5_solution, construct_approximate_NSE_3_5_7_solution, construct_best_approximate_NSE_3_5_7_solution

end

module CoupledGL

using ..Utilities: periodic_time

const RESULTS_DIR = normpath(joinpath(@__DIR__, "..", "results"))

include("analytical_solutions/CoupledGLAnalytical.jl")
using .CoupledGLAnalytical
export CoupledGLAnalytical

include("analytical_solutions/CoupledGLResiduals.jl")
using .CoupledGLResiduals
export CoupledGLResiduals

include("solvers/solve_coupled.jl")
export solve_coupled_gl

include("convergence.jl")
export max_relative_error_over_run, save_convergence_csv

export RESULTS_DIR, periodic_time

end

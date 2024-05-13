module Solvers

using ProgressMeter, LinearAlgebra
using GeneralizedSchrodingerResearch.Utilities: filtration, create_FD_matrix
using GeneralizedSchrodingerResearch.Utilities: integral_1, integral_2

include("general_solver.jl")
export solve

end
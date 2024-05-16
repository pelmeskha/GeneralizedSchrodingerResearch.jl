module Solvers

using ProgressMeter, LinearAlgebra
using GeneralizedSchrodingerResearch.BaseUtilities: filtration, create_FD_matrix
using GeneralizedSchrodingerResearch.BaseUtilities: integral_1, integral_2

include("general_solver.jl")
export solve

end
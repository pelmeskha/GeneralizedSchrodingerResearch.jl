module Solvers

using ProgressMeter
using GeneralizedSchrodingerResearch.Utilities: filtration
using GeneralizedSchrodingerResearch.Utilities: integral_1

include("fourier.jl")
export fourier_solve

end
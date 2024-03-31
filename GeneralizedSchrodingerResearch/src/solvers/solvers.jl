module Solvers

using ProgressMeter
using GeneralizedSchrodingerResearch.Utilities: filtration
using GeneralizedSchrodingerResearch.Utilities: integral_1, integral_2

include("fourier.jl")
export fourier_solve

end
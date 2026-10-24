module GeneralizedSchrodingerResearchGPU

using CUDA
using GLMakie
using LinearAlgebra
using Parameters
using ProgressMeter
using GeneralizedSchrodingerResearch.BaseUtilities: create_FD_matrix
import GeneralizedSchrodingerResearch.BaseUtilities: cuda_integral_1, cuda_integral_2, cuda_filtration
import GeneralizedSchrodingerResearch.GNSE_3_5_7: cuda_solve, live_plot_record

include("cuda/cuda_integrals.jl")
include("cuda/cuda_filtration.jl")
include("cuda/cuda_operations.jl")
include("cuda/cuda_solve.jl")
include("live_plot.jl")

end

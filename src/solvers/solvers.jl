module Solvers

using ProgressMeter, LinearAlgebra, CUDA
using GeneralizedSchrodingerResearch.BaseUtilities: filtration, create_FD_matrix,
    integral_1, integral_2, cuda_integral_1, cuda_integral_2, cuda_filtration

include("solve.jl")
export solve

include("cuda_operations.jl")
export cuda_matrix_vector_multiplication, cuda_simple_tolerance

include("cuda_solve.jl")
export cuda_solve

end
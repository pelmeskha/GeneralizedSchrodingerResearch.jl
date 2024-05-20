module BaseUtilities

using CUDA, Statistics

include("integrals.jl")
export integral_1, integral_2

include("filtration.jl")
export filtration

include("create_FD_matrix.jl")
export create_FD_matrix

include("parameter_transform.jl")
export M0_M1_to_ε2_ε3, ε2_ε3_to_M0_M1

include("cuda_integrals.jl")
export cuda_integral_1, cuda_integral_2

include("cuda_filtration.jl")
export cuda_filtration

include("smooth_vector.jl")
export smooth_vector

end
module BaseUtilities

using CUDA, Statistics, Plots

include("integrals.jl")
export integral_1, integral_2

include("filtration.jl")
export filtration

include("create_FD_matrix.jl")
export create_FD_matrix

include("metrics.jl")
export relative_error_to_amplitude

include("parameter_transform.jl")
export M0_M1_to_ε2_ε3, ε2_ε3_to_M0_M1

include("cuda_integrals.jl")
export cuda_integral_1, cuda_integral_2

include("cuda_filtration.jl")
export cuda_filtration

include("parse_to_vector.jl")
export parse_to_vector

include("plot_utilities.jl")
export set_grid_alpha

include("shift_pulse_to_center.jl")
export shift_pulse_to_center, find_threshold

include("smooth_vector.jl")
export smooth_vector

end
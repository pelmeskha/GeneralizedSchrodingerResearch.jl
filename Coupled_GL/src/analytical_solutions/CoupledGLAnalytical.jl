module CoupledGLAnalytical

include("modeling_args.jl")
include("tanh.jl")
include("sn.jl")
include("bright.jl")
include("vector.jl")

export ModelingArgs
export construct_dependent_args_for_tanh, tanh_amplitude, tanh_solution
export construct_dependent_args_for_sn, sn_parameters, sn_solution, jacobi_sn, complete_elliptic_K
export construct_dependent_args_for_bright, bright_solution
export construct_dependent_args_for_vector, vector_u_solution, vector_v_solution
export coupled_gl_xi

end # module

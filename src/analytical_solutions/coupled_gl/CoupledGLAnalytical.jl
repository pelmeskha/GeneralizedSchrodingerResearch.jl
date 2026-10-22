module CoupledGLAnalytical

include("modeling_args.jl")
include("tanh.jl")
include("sn.jl")
include("bright.jl")

export ModelingArgs
export construct_dependent_args_for_tanh, tanh_amplitude, tanh_solution
export construct_dependent_args_for_sn, sn_parameters, sn_solution, jacobi_sn, complete_elliptic_K
export construct_dependent_args_for_bright, bright_solution
export coupled_gl_xi, field_error

"""
    field_error(U, T, z, true_field)

Compare a numerical field `U` (size `(length(T), length(z))`) with an analytical
`true_field(z, T)` on the same grid.

Returns a named tuple:
- `max_modulus_error`: max over z, T of `| |U| - |true| |`;
- `max_complex_error`: max over z, T of `| U - true |` (also sees phase errors);
- `max_relative_error_percent`: `max_complex_error` as a percentage of `max |true|`
  (`NaN` when the reference is identically zero).
"""
function field_error(U, T, z, true_field)
    ref = [true_field(zi, tj) for tj in T, zi in z]
    scale = maximum(abs, ref)
    max_modulus_error = maximum(abs.(abs.(U) .- abs.(ref)))
    max_complex_error = maximum(abs.(U .- ref))
    max_relative_error_percent = scale > 0 ? 100 * max_complex_error / scale : NaN
    return (
        max_modulus_error=max_modulus_error,
        max_complex_error=max_complex_error,
        max_relative_error_percent=max_relative_error_percent,
    )
end

end # module

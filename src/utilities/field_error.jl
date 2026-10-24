"""
    field_error(U, T, z, true_field)

Compare a numerical field `U` (size `(length(T), length(z))`) with an analytical
reference on the same grid. The reference is either a callable `true_field(z, T)`
or a matrix of size `(length(T), length(z))`.

Returns a named tuple:
- `max_modulus_error`: max over z, T of `| |U| - |true| |`;
- `max_complex_error`: max over z, T of `| U - true |` (also sees phase errors);
- `max_relative_error_percent`: `max_complex_error` as a percentage of `max |true|`
  (`NaN` when the reference is identically zero).
"""
function field_error(U, T, z, true_field::Function)
    return field_error(U, T, z, [true_field(zi, tj) for tj in T, zi in z])
end

function field_error(U, T, z, ref::AbstractMatrix)
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

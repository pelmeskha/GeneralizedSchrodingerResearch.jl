module CoupledGLResiduals

"""
    coupled_gl_residual(
        u_solution, v_solution, z, t;
        Δβ, δ, β2, g, Ω_g, γ,
        hz=1e-3, ht=1e-3
    )

    Evaluate the coupled Ginzburg–Landau PDE residuals at (z, t).

    Solutions must be callable as:
        u_solution(z, t)
        v_solution(z, t)

    Derivatives use fourth-order centered finite differences.
    The functions are evaluated at z ± hz, z ± 2hz,
    and t ± ht, t ± 2ht.

    Returns:
        residual_u, residual_v : complex PDE residuals
        error_u, error_v       : absolute residual magnitudes

    A smooth exact solution should produce residuals near zero,
    up to differentiation and floating-point errors.
"""
function coupled_gl_residual(
    u_solution, v_solution;
    z, t, Δβ, δ, β2, g, Ω_g, γ,
    hz=1e-3,
    ht=1e-3
)
    hz > 0 || throw(ArgumentError("hz must be positive"))
    ht > 0 || throw(ArgumentError("ht must be positive"))
    Ω_g != 0 || throw(ArgumentError("Ω_g must be nonzero"))

    function field_derivatives(field)
        f = field(z, t)

        fz_minus2 = field(z - 2*hz, t)
        fz_minus1 = field(z - hz, t)
        fz_plus1  = field(z + hz, t)
        fz_plus2  = field(z + 2*hz, t)

        ft_minus2 = field(z, t - 2*ht)
        ft_minus1 = field(z, t - ht)
        ft_plus1  = field(z, t + ht)
        ft_plus2  = field(z, t + 2*ht)

        f_z = (
            fz_minus2 - 8*fz_minus1 +
            8*fz_plus1 - fz_plus2
        ) / (12*hz)

        f_t = (
            ft_minus2 - 8*ft_minus1 +
            8*ft_plus1 - ft_plus2
        ) / (12*ht)

        f_tt = (
            -ft_plus2 + 16*ft_plus1 - 30*f +
            16*ft_minus1 - ft_minus2
        ) / (12*ht^2)

        return f, f_z, f_t, f_tt
    end

    u, u_z, u_t, u_tt = field_derivatives(u_solution)
    v, v_z, v_t, v_tt = field_derivatives(v_solution)

    residual_u = (
        u_z
        - 1im * Δβ/2 * u
        + δ * u_t
        + 1im * β2/2 * u_tt
        - g/2 * u
        - g/(2*Ω_g^2) * u_tt
        - 1im * γ * (abs2(u) + (2/3)*abs2(v)) * u
        - 1im * γ/3 * v^2 * conj(u)
    )

    residual_v = (
        v_z
        - 1im * Δβ/2 * v
        + δ * v_t
        + 1im * β2/2 * v_tt
        - g/2 * v
        - g/(2*Ω_g^2) * v_tt
        - 1im * γ * (abs2(v) + (2/3)*abs2(u)) * v
        - 1im * γ/3 * u^2 * conj(v)
    )

    return (
        residual_u=residual_u,
        residual_v=residual_v,
        error_u=abs(residual_u),
        error_v=abs(residual_v)
    )
end
"""
    Evaluate residual magnitudes on a grid.

    Output matrices have shape:
        (length(t_values), length(z_values))
    # Check if equations satisfied
    parameters = (Δβ=Δβ, δ=δ, β2=β2, g=g, Ω_g=Ω_g, γ=γ)

    check = coupled_gl_residual_grid(
        true_u,
        true_v;
        z_values = range(0.0, 40.0, length=101),
        t_values = range(Tspan[1], Tspan[2], length=256),
        parameters...,
        hz=1e-4,
        ht=1e-4
    )

    println("True solution check:")
    println("Maximum |R_u|: ", check.max_error_u)
    println("Maximum |R_v|: ", check.max_error_v)
"""
function coupled_gl_residual_grid(
    u_solution, v_solution;
    z_values, t_values,
    kwargs...
)
    error_u = Matrix{Float64}(
        undef, length(t_values), length(z_values)
    )
    error_v = similar(error_u)

    for (i, zi) in enumerate(z_values)
        for (j, tj) in enumerate(t_values)
            result = coupled_gl_residual(
                u_solution, v_solution;
                z=zi, t=tj,
                kwargs...
            )

            error_u[j, i] = result.error_u
            error_v[j, i] = result.error_v
        end
    end

    return (
        error_u=error_u,
        error_v=error_v,
        max_error_u=maximum(error_u),
        max_error_v=maximum(error_v)
    )
end

export coupled_gl_residual, coupled_gl_residual_grid

end # module

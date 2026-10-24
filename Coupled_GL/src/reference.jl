function sn_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    k = 12.6
    c0 = 1.4
    mu = 0.97
    I1 = 0.1
    ω = 0
    w0 = 0

    β2 = 0.027168436376788473
    Ω_g = 66.08752874553595
    Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω) / (2*c0)
    γ  = mu^4 * c0^2 * β2 / (32*I1)

    args = ModelingArgs(
        (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
        construct_dependent_args_for_sn,
    )

    u0(T) = sn_solution(0.0, T, args)
    v0(T) = 0

    true_u(z, T) = sn_solution(z, T, args)
    true_v(z, T) = 0

    true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, c0, w0))
    true_v_periodic(z, T) = true_v(z, periodic_time(z, T, Tspan, c0, w0))

    return (; tag="sn", Tspan, zspan, args, u0, v0, true_u, true_v, true_u_periodic, true_v_periodic)
end

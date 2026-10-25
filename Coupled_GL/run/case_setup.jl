using GeneralizedSchrodingerResearch.Utilities: periodic_time
using GeneralizedSchrodingerResearch.CoupledGL.CoupledGLAnalytical: ModelingArgs, construct_dependent_args_for_sn, sn_solution, construct_dependent_args_for_tanh, tanh_solution, construct_dependent_args_for_bright, bright_solution, construct_dependent_args_for_vector, vector_u_solution, vector_v_solution

function build_reference(tag, Tspan, zspan, args, solution)
    u0(T) = solution(0.0, T, args)
    v0(T) = 0

    true_u(z, T) = solution(z, T, args)
    true_v(z, T) = 0

    true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, args.c0, args.w0))
    true_v_periodic(z, T) = true_v(z, periodic_time(z, T, Tspan, args.c0, args.w0))

    return (; tag, Tspan, zspan, args, u0, v0, true_u, true_v, true_u_periodic, true_v_periodic)
end

function reference()
    return sn_reference()
end

function tanh_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13, ω=0, β2=0.027, Ω_g=66, w0=0, Tspan=Tspan, mu=0.97, I1=0.1),
        construct_dependent_args_for_tanh,
    )

    return build_reference("tanh", Tspan, zspan, args, tanh_solution)
end

function sn_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    k = 13
    c0 = 1.4
    mu = 0.97
    I1 = 0.1
    ω = 0
    w0 = 0

    β2 = 0.027
    Ω_g = 66.0
    Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω) / (2*c0)
    γ  = mu^4 * c0^2 * β2 / (32*I1)

    args = ModelingArgs(
        (c0=c0, k=k, ω=ω, β2=β2, Ω_g=Ω_g, w0=w0, Tspan=Tspan, Δβ=Δβ, γ=γ),
        construct_dependent_args_for_sn,
    )

    return build_reference("sn", Tspan, zspan, args, sn_solution)
end

function bright_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13, ω=5, β2=0.025, Ω_g=10, w0=0, Tspan=Tspan, Δβ=25, γ=-0.014),
        construct_dependent_args_for_bright,
    )

    return build_reference("bright", Tspan, zspan, args, bright_solution)
end

function build_vector_reference(tag, Tspan, zspan, args)
    u0(T) = vector_u_solution(0.0, T, args)
    v0(T) = vector_v_solution(0.0, T, args)

    true_u(z, T) = vector_u_solution(z, T, args)
    true_v(z, T) = vector_v_solution(z, T, args)

    true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, args.c0, args.w0))
    true_v_periodic(z, T) = true_v(z, periodic_time(z, T, Tspan, args.c0, args.w0))

    return (; tag, Tspan, zspan, args, u0, v0, true_u, true_v, true_u_periodic, true_v_periodic)
end

function vector_linear_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13.52, ω=1.5, β2=0.027, Ω_g=66, w0=0, Tspan=Tspan, Δβ=25, γ=-0.07, ψ0=0.0, θ=π/8),
        construct_dependent_args_for_vector,
    )

    return build_vector_reference("vector_linear", Tspan, zspan, args)
end

function vector_circular_reference()
    Tspan = (-10.0, 10.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13.52, ω=1.5, β2=0.027, Ω_g=66, w0=0, Tspan=Tspan, Δβ=25, γ=-0.07, ψ0=π/2, θ=0.0),
        construct_dependent_args_for_vector,
    )

    return build_vector_reference("vector_circular", Tspan, zspan, args)
end

# Non-periodic domain for the finite-difference solver. The pulse travels from T = 0 at z = 0
# to T ≈ 36 at z = 50, so the window must contain it with margin for the field to be negligible at the ends.
# Kink solution u = A tanh(μ(z - c0 T)/2) on the non-periodic window; the field tends to ±A at the ends,
# so it needs the Neumann boundary condition rather than the Dirichlet one.
function tanh_open_reference()
    Tspan = (-25.0, 60.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13, ω=0, β2=0.027, Ω_g=66, w0=0, Tspan=Tspan, mu=0.97, I1=0.1),
        construct_dependent_args_for_tanh,
    )

    true_u(z, T) = tanh_solution(z, T, args)
    true_v(z, T) = 0
    return (; tag="tanh_open", Tspan, zspan, args,
        u0=T -> true_u(0.0, T), v0=T -> 0.0, true_u, true_v)
end

function vector_circular_open_reference()
    Tspan = (-25.0, 60.0)
    zspan = (0.0, 50.0)

    args = ModelingArgs(
        (c0=1.4, k=13.52, ω=1.5, β2=0.027, Ω_g=66, w0=0, Tspan=Tspan, Δβ=25, γ=-0.07, ψ0=π/2, θ=0.0),
        construct_dependent_args_for_vector,
    )

    true_u(z, T) = vector_u_solution(z, T, args)
    true_v(z, T) = vector_v_solution(z, T, args)
    return (; tag="vector_circular_open", Tspan, zspan, args,
        u0=T -> true_u(0.0, T), v0=T -> true_v(0.0, T), true_u, true_v)
end

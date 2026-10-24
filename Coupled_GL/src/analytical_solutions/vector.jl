"""
    construct_dependent_args_for_vector(indep)

Required independent keys: `c0, k, ω, β2, Δβ, γ, w0, ψ0, θ`.
Derived: `eta, A, B` (η² = −ξ₁; with X = β₂c₀²ξ₁/γ the linear family cos 2ψ₀ > 0 has A = √X cos θ, B = √X sin θ,
the circular family cos 2ψ₀ < 0 has A = B = √(3X/4)); `δ, g` from the common restrictions.
`θ` is used only by the linear family. Requires ξ₁ < 0, β₂γ < 0 and sin 2ψ₀ = 0, otherwise throws `ArgumentError`.
"""
function construct_dependent_args_for_vector(indep)
    (; c0, k, ω, β2, Δβ, γ, ψ0, θ) = indep
    xi = coupled_gl_xi(; Δβ, β2, γ, c0, k, ω)
    xi.xi1 < 0 || throw(ArgumentError("vector soliton requires ξ1 < 0, got ξ1 = $(xi.xi1)"))
    β2 * γ < 0 || throw(ArgumentError("vector soliton requires β2*γ < 0, got β2*γ = $(β2*γ)"))
    abs(sin(2ψ0)) < 1e-12 || throw(ArgumentError("vector soliton requires sin(2ψ0) = 0, got ψ0 = $ψ0"))
    X = β2 * c0^2 * xi.xi1 / γ
    if cos(2ψ0) > 0
        A = sqrt(X) * cos(θ)
        B = sqrt(X) * sin(θ)
    else
        A = sqrt(3X / 4)
        B = A
    end
    return merge(common_dependent_args(indep), (eta=sqrt(-xi.xi1), A=A, B=B))
end

function vector_u_solution(z, t, args::ModelingArgs)
    return args.A * sech(args.eta * (z - args.c0*t - args.w0)) * exp(1im * (args.k*z - args.ω*t))
end

function vector_v_solution(z, t, args::ModelingArgs)
    return args.B * sech(args.eta * (z - args.c0*t - args.w0)) * exp(1im * (args.k*z - args.ω*t + args.ψ0))
end

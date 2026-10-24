"""
    construct_dependent_args_for_bright(indep)

Required independent keys: `c0, k, ω, β2, Ω_g, w0, Δβ, γ`.
Derived: `eta, A` (η² = −ξ₁, A² = −2ξ₁/ξ₂); `δ, g` from the common restrictions.
Requires ξ₁ < 0 and ξ₂ > 0, otherwise throws `ArgumentError`.
"""
function construct_dependent_args_for_bright(indep)
    (; c0, k, ω, β2, Δβ, γ) = indep
    xi = coupled_gl_xi(; Δβ, β2, γ, c0, k, ω)
    xi.xi1 < 0 || throw(ArgumentError("bright soliton requires ξ1 < 0, got ξ1 = $(xi.xi1)"))
    xi.xi2 > 0 || throw(ArgumentError("bright soliton requires ξ2 > 0, got ξ2 = $(xi.xi2)"))
    println("ξ1 = $(xi.xi1), ξ2 = $(xi.xi2)")
    return merge(common_dependent_args(indep), (eta=sqrt(-xi.xi1), A=sqrt(-2*xi.xi1/xi.xi2)))
end

function bright_solution(z, t, args::ModelingArgs)
    return args.A * sech(args.eta * (z - args.c0*t - args.w0)) * exp(1im * (args.k*z - args.ω*t))
end

tanh_amplitude(I1, mu) = 2 * sqrt(2 * I1) / abs(mu)

"""
    construct_dependent_args_for_tanh(indep)

Required independent keys: `c0, k, ω, β2, Ω_g, w0, mu, I1`.
Derived: `Δβ, γ` from ξ₁ = μ²/2 and ξ₂ = −μ⁴/(16 I₁); `A`; `δ, g` from the common restrictions.
"""
function construct_dependent_args_for_tanh(indep)
    (; c0, k, ω, β2, mu, I1) = indep
    γ = mu^4 * c0^2 * β2 / (32 * I1)
    Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω) / (2*c0)
    return merge(common_dependent_args(indep), (Δβ=Δβ, γ=γ, A=tanh_amplitude(I1, mu)))
end

function tanh_solution(z, t, args::ModelingArgs)
    return args.A * tanh(args.mu * (z - args.c0*t - args.w0) / 2) * exp(1im * (args.k*z - args.ω*t))
end

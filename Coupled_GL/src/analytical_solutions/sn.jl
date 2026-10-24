function agm_sequences(kp)
    a = [1.0]
    b = [kp]
    c = [sqrt(1 - kp^2)]
    for _ in 1:60
        abs(c[end]) < 1e-17 && break
        push!(a, (a[end] + b[end]) / 2)
        push!(c, (a[end-1] - b[end]) / 2)
        push!(b, sqrt(a[end-1] * b[end]))
    end
    return a, b, c
end

complete_elliptic_K(kp) = π / (2 * agm_sequences(kp)[1][end])

jacobi_sn(u, kp::Real) = jacobi_sn(u, agm_sequences(kp))

function jacobi_sn(u, agm::Tuple)
    a, _, c = agm
    N = length(a) - 1
    phi = 2.0^N * a[end] * u
    for n in N:-1:1
        phi = (phi + asin(clamp(c[n+1] / a[n+1] * sin(phi), -1.0, 1.0))) / 2
    end
    return sin(phi)
end

function sn_parameters(xi1, xi2, period)
    amplitude_and_kappa(kp) = (
        sqrt(2 * (1 - kp^2) * xi1 / ((2 - kp^2) * abs(xi2))),
        sqrt(xi1 / (2 - kp^2)),
    )
    lo, hi = log(1e-300), 0.0
    for _ in 1:200
        mid = (lo + hi) / 2
        kp = exp(mid)
        _, kappa = amplitude_and_kappa(kp)
        if 4 * complete_elliptic_K(kp) / kappa > period
            lo = mid
        else
            hi = mid
        end
    end
    kp = exp((lo + hi) / 2)
    amplitude, kappa = amplitude_and_kappa(kp)
    return (A = amplitude, kappa = kappa, kp = kp)
end

"""
    construct_dependent_args_for_sn(indep)

Required independent keys: `c0, k, ω, β2, Ω_g, w0, Δβ, γ, Tspan`.
Derived: `A, kappa, kp` (from ξ₁, ξ₂ and the period c₀·(Tspan[2] − Tspan[1])); `agm` (AGM sequences of `kp`,
computed once for fast evaluation); `δ, g` from the common restrictions.
"""
function construct_dependent_args_for_sn(indep)
    (; c0, k, ω, β2, Δβ, γ, Tspan) = indep
    xi = coupled_gl_xi(; Δβ, β2, γ, c0, k, ω)
    pars = sn_parameters(xi.xi1, xi.xi2, c0 * (Tspan[2] - Tspan[1]))
    return merge(common_dependent_args(indep), pars, (agm = agm_sequences(pars.kp),))
end

function sn_solution(z, t, args::ModelingArgs)
    return args.A * jacobi_sn(args.kappa * (z - args.c0*t - args.w0), args.agm) * exp(1im * (args.k*z - args.ω*t))
end

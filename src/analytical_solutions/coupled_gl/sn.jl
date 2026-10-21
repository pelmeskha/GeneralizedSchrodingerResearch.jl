coupled_gl_xi(; Δβ, β2, γ, c0, k, ω) = (
    xi1 = (-Δβ - 2*ω/c0 + β2*ω^2 + 2*k) / (β2*c0^2),
    xi2 = -2*γ / (β2*c0^2),
)

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

function jacobi_sn(u, kp)
    a, _, c = agm_sequences(kp)
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

function sn_solution(z, t, pars; c0, w0, k, ω)
    return pars.A * jacobi_sn(pars.kappa * (z - c0*t - w0), pars.kp) * exp(1im * (k*z - ω*t))
end

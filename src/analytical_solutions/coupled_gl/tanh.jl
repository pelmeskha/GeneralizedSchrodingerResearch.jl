tanh_amplitude(I1, mu) = 2 * sqrt(2 * I1) / abs(mu)

function tanh_solution(z, t; I1, mu, c0, w0, k, ω)
    A = tanh_amplitude(I1, mu)
    return A * tanh(mu * (z - c0*t - w0) / 2) * exp(1im * (k*z - ω*t))
end

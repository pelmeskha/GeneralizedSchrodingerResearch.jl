function NSE_5_soliton(x, t, k, ω, ε_2, theta_0, z_0)
    μ = 4*(ω-k^2)
    μ > 0.0 || throw(ArgumentError("μ ≤ 0. Check k and ω."))
    power = sqrt(μ) * (x - 2k * t - z_0)
    ν = 4/3 * ε_2
    Complex(
        sqrt(
            (4 * μ * exp(power)) /
            (1 + 4 * exp(power) + (4 + 4 * μ * ν) * exp(2 * power))
        )
    ) *
    exp(
        1im * (k * x - theta_0 - ω*t)
    )
end
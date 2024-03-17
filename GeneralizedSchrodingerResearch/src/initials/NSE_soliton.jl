function NSE_soliton(x, t, k, ω, theta_0, z_0)
    Complex(
        4 * (k^2 - ω) /
        (
            2 * (
                k^2 - ω
            ) * exp(
                -(x - z_0 - 2k * t) * sqrt(k^2 - ω)
            ) + exp(
                (x - z_0 - 2k * t) * sqrt(k^2 - ω)
            )
        )
    ) *
    exp(
        1im * (k * x - theta_0 - ω*t)
    )
end
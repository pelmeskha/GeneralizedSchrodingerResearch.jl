function NSE_soliton(x, b_1, k, omega, theta_0, z_0)
    Complex.(
        4 * (k^2 - omega) ./
        (
            2 * b_1 * (
                k^2 - omega
            ) * exp.(
                -(x .- z_0) * sqrt(k^2 - omega)
            ) + exp.(
                (x .- z_0) * sqrt(k^2 - omega)
            )
        )
    ) .*
    exp.(
        1im * (k * x .- theta_0)
    )
end
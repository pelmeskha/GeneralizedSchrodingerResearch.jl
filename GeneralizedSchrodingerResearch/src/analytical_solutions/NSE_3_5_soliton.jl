function NSE_3_5_soliton(
    x,
    t,
    k::Real,
    ω::Real,
    ε_2::Real,
    theta_0::Real,
    z_0::Real; 
    cycle::Bool=false,
    L::Real=0.0,
)
    μ = 4*(ω-k^2)
    μ > 0.0 || throw(ArgumentError("μ ≤ 0. Check k and ω."))
    c=2*k
    if cycle
        if t>(L/2+x)/c
            t-=L/c*floor(1/2 + (c*t -x)/L)
        end
    end
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
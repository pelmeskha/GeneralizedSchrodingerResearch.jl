function NSE_3_soliton(
    x,
    t,
    k::Real,
    ω::Real,
    theta_0::Real,
    z_0::Real;
    cycle::Bool=false,
    L::Real=0.0,
)
    μ = (ω-k^2) # Original paper: k^2 - ω
    μ > 0.0 || throw(ArgumentError("μ ≤ 0. Check k and ω."))
    c=2*k
    if cycle
        if t>(L/2+x)/c
            t-=L/c*floor(1/2 + (c*t -x)/L)
        end
    end
    Complex(
        4 * (μ) / 
        (
            2 * (
                μ
            ) * exp(
                -(x - z_0 - 2k * t) * sqrt(μ)
            ) + exp(
                (x - z_0 - 2k * t) * sqrt(μ)
            )
        )
    ) *
    exp(
        1im * (k * x - theta_0 - ω*t)
    )
end
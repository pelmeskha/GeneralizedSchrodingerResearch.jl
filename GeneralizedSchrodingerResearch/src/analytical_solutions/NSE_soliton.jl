function NSE_soliton(
    x,
    t,
    k,
    ω,
    theta_0,
    z_0;
    cycle::Bool=false,
    L::Float64=0.0,
    c::Float64=0.0,
)
    μ = (ω-k^2) # Original paper: k^2 - ω
    μ > 0.0 || throw(ArgumentError("μ ≤ 0. Check k and ω."))

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
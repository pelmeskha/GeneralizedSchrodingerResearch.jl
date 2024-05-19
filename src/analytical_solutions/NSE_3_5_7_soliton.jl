function evaluate_μ(M₀, M₁)
    μ_a = sqrt(M₁ / (M₀ * (M₁ + 6 * M₀)))
    μ_b = -μ_a #TODO проверить влияние
    return μ_a
end
function ξ_edges(μ, M₀, M₁, ε_left, ε_right, ξ₀)
    sqrt_expression = sqrt(4 * M₀ * M₁ + M₁^2)
    return (
        log((sqrt_expression - 2 * M₀ - M₁) / (2 * M₀)) / μ + ξ₀ + ε_left,
        log((-sqrt_expression - 2 * M₀ - M₁) / (2 * M₀)) / μ + ξ₀ - ε_right
    )
end
function hyperbolic_space(a::Real, b::Real, N::Int; density::Real=1.0)
    N ≥ 2 || error("You must specify at least 2 points")
    density > 0 || error("Parameter density must be positive")
    linear_space = range(-1.0, stop=1.0, length=N)
    hyperbolic_space = tanh.(linear_space * 2 * density) / tanh(2 * density)
    return a .+ (b - a) .* (hyperbolic_space .+ 1) / 2
end
function evaluate_z(ξ_vector, μ, M₀, M₁, z₀, ξ₀)
    sqrt_expression = sqrt(4 * M₀ * M₁ + M₁^2)
    _atanh_expression = (2 * M₀ * exp.(μ * (ξ_vector .- ξ₀)) .+ 2 * M₀ .+ M₁) / sqrt_expression
    v_clipped = [x > 1 ? 1 : x < -1 ? -1 : x for x in _atanh_expression]
    atanh_expression = unique(v_clipped)
    return z₀ .+ ξ_vector / M₀ + (2 * M₁) / (μ * M₀ * sqrt_expression) *
        atanh.(atanh_expression)
end
function evaluate_y(ξ_vector, μ, M₀, M₁, ξ₀)
    exp_expression = (1 .+ exp.(μ * (ξ_vector .- ξ₀)))
    
    complex_y = sqrt.(Complex.(M₀ .+ M₁ ./ exp_expression - M₁ ./ (exp_expression .^ 2)))
    all(imag(complex_y).==0) || @warn "При вычислении y(z) обнаружены и отброшены комплесные числа."
    complex_y[imag(complex_y).!=0.0].=Complex(0.0)
    return Float64[real(x) for x in complex_y]
end
function precompile_NSE_3_5_7_soliton(
    ε₂::Real,
    ε₃::Real,
    z₀::Real,
    ξ₀::Real,
    L::Real,
)
    (M₀, M₁) = ε2_ε3_to_M0_M1(ε₂, ε₃)
    println("M₀=",round(M₀,digits=4)," M₁=",round(M₁,digits=4))
    μ = evaluate_μ(M₀, M₁)

    iters, iters_limit, success, ε_left, ε_right = 0, 50, false, 0, 0
    while ~success && iters<iters_limit
        iters+=1
        (ξ_left, ξ_right) = ξ_edges(μ, M₀, M₁, ε_left, ε_right, ξ₀)
        global ξ_vector = hyperbolic_space(ξ_left, ξ_right, 10^6; density=1)
        global z_vector = evaluate_z(ξ_vector, μ, M₀, M₁, z₀, ξ₀)
        if z_vector[1]==-Inf
            ε_left+=eps()
        end
        if z_vector[end]==Inf
            ε_right+=eps()
        end
        success = ~any(isinf, z_vector)
        iters==iters_limit && error("solution precompilation failed: MaxIters")
    end
    maximum(z_vector) > L/2 || @warn "solution precompilation: BoundsError on right, extrapolation will be used"
    minimum(z_vector) < -L/2 || @warn "solution precompilation: BoundsError on left, extrapolation will be used"

    y_vector = evaluate_y(ξ_vector, μ, M₀, M₁, ξ₀)
    println(
        "interpolation defined for [",
        round(minimum(z_vector),digits=4),
        ", ",
        round(maximum(z_vector),digits=4),
        "] with center z=",
        z_vector[argmax(y_vector)]
    )
    interpolator = extrapolate(interpolate((z_vector,), y_vector, Gridded(Linear())), Line())
    return interpolator
end
function reduce_negative_values(x::Real)::Real
    value = interpolator(x)
    return value < 0 ? 0 : value
end
function reduce_negative_values(xs::AbstractVector{<:Real})::Vector{Real}
    return [reduce_negative_values(x) for x in xs]
end
function NSE_3_5_7_soliton(
    x,
    t,
    k::Real,
    ω::Real,
    θ₀::Real,
    z₀::Real,
    precompiled_data;
    # Cycling parameters
    cycle::Bool=false,
    L::Real=0.0,
)
    y = precompiled_data
    c=2*k
    if cycle
        if t > (L / 2 + x) / c
            t -= L / c * floor(1 / 2 + (c * t - x) / L)
        end
    end
    z = x - 2k * t - z₀
    
    reduce_negative_values(y.(z)) * exp(
        1im * (k * x - θ₀ - ω * t)
    )
end
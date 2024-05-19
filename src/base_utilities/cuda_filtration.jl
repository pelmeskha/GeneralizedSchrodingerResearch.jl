function cuda_filtration(
    U::Union{CuArray{ComplexF64, 1}, CuArray{Float64, 1}},
    h::Float64,
    factor::Float64,
    l_nominal::Float64,
)
    delta = trunc(Int, l_nominal / h / 2)
    N = size(U, 1)
    N != 1 || throw(BoundsError("solution N = 1. Consider transpose(U)."))

    i_center = CUDA.argmax(abs.(U))
    i_left = i_center - delta
    i_right = i_center + delta

    i_left = i_left < 1 ? i_left + N : i_left
    i_right = i_right > N ? i_right - N : i_right

    I1 = cuda_integral_1(U, h)
    I2 = cuda_integral_2(U, h)

    if i_right < i_left
        U[i_right:i_left] .= U[i_right:i_left] ./ factor
    else
        U[1:i_left] .= U[1:i_left] ./ factor
        U[i_right:end] .= U[i_right:end] ./ factor
    end

    return U, (I1 - cuda_integral_1(U, h), I2 - cuda_integral_2(U, h))
end
function cuda_matrix_vector_multiplication(M::CuArray{ComplexF64, 2}, x::CuArray{ComplexF64, 1})
    m, _ = size(M)
    y = CUDA.fill(0.0 + 0.0im, m)

    @cuda threads=256 blocks=ceil(Int, m/256) cuda_matrix_vector_multiplication_kernel!(M, x, y)

    return y
end
function cuda_matrix_vector_multiplication_kernel!(M, x, y)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(y)
        sum = 0.0 + 0.0im
        for j in axes(M, 2)
            sum += M[i, j] * x[j]
        end
        y[i] = sum
    end
    return
end
function cuda_calculate_V_kernel!(U, result, tau, ε_2, ε_3)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(U)
        abs_U = abs(U[i])
        phase_factor = exp(1im * tau * (abs_U^2 + ε_2 * abs_U^4 + ε_3 * abs_U^6))
        result[i] = phase_factor * U[i]
    end
    return
end
function cuda_calculate_V(U::CuArray{ComplexF64, 1}, tau::Float64, ε_2::Float64, ε_3::Float64)
    result = CUDA.fill(0.0 + 0.0im, length(U))
    @cuda threads=256 blocks=ceil(Int, length(U)/256) cuda_calculate_V_kernel!(U, result, tau, ε_2, ε_3)
    return result
end
function cuda_simple_tolerance(U::CuArray{ComplexF64, 1}, threshold::Real)
    abs_U = abs.(U)
    max_abs_value = CUDA.reduce(max, abs_U)
    result = (max_abs_value - threshold) / max_abs_value * 100.0
    return result
end
function cuda_maximum(U::CuArray{ComplexF64, 1})
    abs_U = abs.(U)
    return CUDA.reduce(max, abs_U)
end
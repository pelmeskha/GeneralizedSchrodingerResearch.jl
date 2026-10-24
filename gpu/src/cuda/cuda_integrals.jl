function derivative_kernel!(U, h, dU)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(U)
        if i == 1
            dU[i] = (U[2] - U[end]) / (2.0 * h)
        elseif i == length(U)
            dU[i] = (U[1] - U[end-1]) / (2.0 * h)
        else
            dU[i] = (U[i+1] - U[i-1]) / (2.0 * h)
        end
    end
    return
end
function integral_1_kernel!(U, h, result)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(U)
        result[i] = h * abs(U[i])^2
    end
    return
end
function integral_2_kernel!(U, h, dU_dx, result)
    i = threadIdx().x + (blockIdx().x - 1) * blockDim().x
    if i <= length(U)
        dU_dx[i] = dU_dx[i] / h
        result[i] = real(-1im * h * (dU_dx[i] * conj(U[i]) - U[i] * conj(dU_dx[i])))
    end
    return
end
function cuda_integral_1(U::CuArray{ComplexF64}, h::Float64)
    result = CUDA.fill(0.0, length(U))

    @cuda threads=256 blocks=ceil(Int, length(U)/256) integral_1_kernel!(U, h, result)

    return sum(result)
end
function cuda_integral_2(U::CuArray{ComplexF64}, h::Float64)
    dU = CUDA.fill(0.0 + 0.0im, length(U))
    dU_dx = CUDA.fill(0.0 + 0.0im, length(U))
    result = CUDA.fill(0.0, length(U))

    @cuda threads=256 blocks=ceil(Int, length(U)/256) derivative_kernel!(U, h, dU)
    @cuda threads=256 blocks=ceil(Int, length(U)/256) integral_2_kernel!(U, h, dU, result)

    return sum(result)
end
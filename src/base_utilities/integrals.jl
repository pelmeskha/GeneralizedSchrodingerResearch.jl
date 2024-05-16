function integral_1(
    U::Union{Vector{ComplexF64}, Vector{Float64}},
    h::Float64,
)
    return h*sum((abs.(U)).^2)
end

#TODO: Optimize
function integral_2(
    U::Union{Vector{ComplexF64}, Vector{Float64}},
    h::Float64,
)
    dU = [(U[2]-U[end]), (U[3:end]-U[1:end-2])..., (U[1]-U[end-1])]./(2.0 * h)
    dU_dx=dU./(ones(length(U))*h)

    return real(
        -1im*h*sum(
            dU_dx .* conj(U) .- U .* conj(dU_dx)
        )
    )
end
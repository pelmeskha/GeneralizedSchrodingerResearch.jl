function integral_1(U,h)
    return h*sum((abs.(U)).^2)
end

function integral_2(U,h)

    append!(U,U[1])
    dU = diff(U)
    pop!(U)

    dU_dx=dU./(ones(length(U))*h)

    return real(
        -1im*h*sum(
            dU_dx .* conj(U) .- U .* conj(dU_dx)
        )
    )
end
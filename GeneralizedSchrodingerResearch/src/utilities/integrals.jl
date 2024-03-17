function integral_1(U,h)
    return h*sum((abs.(U)).^2)
end
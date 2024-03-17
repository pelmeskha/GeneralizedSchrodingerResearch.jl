function filtration(
    U,
    x,
    factor,
    l_nominal,
)
    h = x[2]-x[1]
    delta = trunc(Int, l_nominal / h / 2) # полупротяжённость в ячейках сетки
    N = size(U)[1]
    N != 1 || throw(BoundsError("solution N = 1. Consider transpose(U)."))

    i_center = argmax(abs.(U))
    i_left = i_center - delta
    i_right = i_center + delta

    i_left = i_left < 1 ? i_left + N : i_left
    i_right = i_right > N ? i_right - N : i_right

    I1 = integral_1(U, h)

    if i_right < i_left
        U[i_right:i_left]/=factor
    else
        U[1:i_left]/=factor
        U[i_right:end]/=factor
    end

    return U, I1 - integral_1(U, h)
end
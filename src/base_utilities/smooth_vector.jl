"""
    Функция для сглаживания вектора y, заменяя каждые n значений на их среднее.

    :param y: Входной вектор y.
    :param n: Размер окна для вычисления среднего.
    :return: Сглаженный вектор y.
"""
function smooth_vector(y::Vector{<:Real}, n::Int)
    n < length(y) || throw(ArgumentError("Размер окна n не может быть больше длины вектора y."))
    smoothed_y = copy(y)
    for i in 1:(length(y) - n + 1)
        smoothed_y[i:i+n-1] .= mean(y[i:i+n-1])
    end
    return smoothed_y
end
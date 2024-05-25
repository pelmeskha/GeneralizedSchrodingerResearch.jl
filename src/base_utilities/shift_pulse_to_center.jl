function find_threshold(v::Vector{Float64}, N::Int)
    sorted_v = sort(v, rev=true)
    N <= length(v) || throw(ArgumentError("N = $N не может быть больше длины вектора $(length(v))"))
    threshold_value = sorted_v[N]
    max_value = maximum(v)
    return max_value - threshold_value
end
function shift_pulse_to_center(x::Vector{Float64}, y::Vector{Float64}; y_threshold::Float64 = 0.0)
    length(x) == length(y) || throw(ArgumentError("Vectors x and y must have the same length"))
    max_index = argmax(y)
    n_circ_shift = -max_index + cld(length(x),2) + 1
    y_centered = circshift(y, n_circ_shift)
    filter_indices = findall(y_centered .≥ y_threshold)
    x_filtered = x[filter_indices]
    y_filtered = y_centered[filter_indices]
    return x_filtered, y_filtered, n_circ_shift
end
#TODO: вынести общий код в отдельную функцию
function construct_approximate_NSE_3_5_solution(
    x_grid::Vector{<:Real},
    h::Real,
    U::Union{Vector{ComplexF64},Vector{<:Real}},
    ε₂::Real,
    L::Real;
    debug_flag=false,
)
    abs_U = abs.(U)
    debug_flag && println("Индекс максимума исходного решения: $(argmax(abs_U))")

    # Точное определение максимума (снижение сеточной погрешности B-Сплайновой интерполяцией)
    N_interpolation_points = 12
    bottom_gap = find_threshold(abs_U, N_interpolation_points)
    x_reduced, y_reduced, n_circ_shift = shift_pulse_to_center(x_grid, abs_U; y_threshold=maximum(abs_U) - bottom_gap)

    debug_flag && println("Для центрирования импульс будет сдвинут на $n_circ_shift узлов")

    x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
    spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
    minimization_function = xi -> -spline_interpolator(xi)
    result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
    x_max, y_max = result.minimizer, -result.minimum

    debug_flag && println("Солитон передвинут. Считается, что его интерполированный максимум находится в точке x=$x_max")
    debug_flag && println("Интерполированный максимум = $y_max")

    k = 0.0
    possible_μ = 2 / 3 * (2 * ε₂ * y_max^4 + 3 * y_max^2)
    _possible_ω = possible_μ / 4
    possible_ω = _possible_ω + k^2

    shifted_grid_solution_3_5 = abs.(NSE_3_5_soliton.(x_grid, 0.0, k, possible_ω, ε₂, 0.0, 0.0))

    # x-coordinate correctrion
    bottom_gap_shifted = find_threshold(shifted_grid_solution_3_5, N_interpolation_points)
    x_reduced, y_reduced, _ = shift_pulse_to_center(
        x_grid,
        shifted_grid_solution_3_5;
        y_threshold=maximum(shifted_grid_solution_3_5) - bottom_gap_shifted
    )
    x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
    spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
    minimization_function = xi -> -spline_interpolator(xi)
    result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
    x_max_shifted, _ = result.minimizer, -result.minimum
    _z_0 = -(x_max_shifted - x_max)

    debug_flag && println("Солитон построен, но теперь его интерполированный максимум \
        находится в точке x=$x_max_shifted. Переносим на $_z_0")

    possible_solution_3_5 = x -> NSE_3_5_soliton.(x, 0.0, k, possible_ω, ε₂, 0.0, _z_0)
    deshifted_grid_solution_3_5 = abs.(possible_solution_3_5.(x_grid))
    if debug_flag
        # Проверяем что передвинули куда надо
        bottom_gap_deshifted = find_threshold(deshifted_grid_solution_3_5, N_interpolation_points)
        x_reduced, y_reduced, _ = shift_pulse_to_center(
            x_grid,
            deshifted_grid_solution_3_5;
            y_threshold=maximum(deshifted_grid_solution_3_5) - bottom_gap_deshifted
        )
        x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
        spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
        minimization_function = xi -> -spline_interpolator(xi)
        result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
        x_max_deshifted, _ = result.minimizer, -result.minimum
        println("Солитон перестроен, и теперь его интерполированный максимум находится в точке x=$x_max_deshifted")
        println("Для верности, индекс максимума: $(argmax(deshifted_grid_solution_3_5))")
        println("Добавка для него $(argmax(abs_U)-argmax(deshifted_grid_solution_3_5))")
    end
    return circshift(possible_solution_3_5.(x_grid), argmax(abs_U) - argmax(deshifted_grid_solution_3_5))
end
function construct_approximate_NSE_3_5_7_solution(
    x_grid::Vector{<:Real},
    h::Real,
    abs_U::Vector{<:Real},
    ε₂::Real,
    L::Real;
    use_M_values=false,
    M0::Real=0.0,
    use_M1_value=false,
    M1::Real=0.0,
    debug_flag=false,
    warn_ignore=false,
)

    debug_flag && println("Индекс максимума исходного решения: $(argmax(abs_U))")
    _, y_shifted, n_circ_shift = shift_pulse_to_center(x_grid, abs_U)
    debug_flag && println("Для центрирования импульс будет сдвинут на $n_circ_shift узлов")

    # Точное определение максимума (снижение сеточной погрешности B-Сплайновой интерполяцией)
    N_interpolation_points = 12
    bottom_gap = find_threshold(abs_U, N_interpolation_points)
    x_reduced, y_reduced, _ = shift_pulse_to_center(x_grid, abs_U; y_threshold=maximum(abs_U) - bottom_gap)
    x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
    spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
    minimization_function = xi -> -spline_interpolator(xi)
    result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
    x_max, y_max = result.minimizer, -result.minimum
    debug_flag && println("Солитон передвинут. Считается, что его интерполированный максимум находится в точке x=$x_max")
    debug_flag && println("Интерполированный максимум = $y_max")
    if use_M_values
        if ~use_M1_value
            M1 = 4 * (y_max^2 - M0)
        end
        debug_flag && println("M1=$M1")
        _ε₂, _ε₃ = 0.0, 0.0
    else
        D = 2 * sqrt(9 + 4 * ε₂^2 * y_max^4 - 6 * y_max^2 * ε₂)
        b = -6 - 4 * ε₂ * y_max^2
        M0 = (b - D) / (4 * ε₂)
        M1 = 4 * (y_max^2 - M0)
        _ε₂, _ε₃ = M0_M1_to_ε2_ε3(M0, M1)
        debug_flag && println("Параметры подобраны для ε₂ = ", _ε₂, ", ε₃ = ", _ε₃)
    end

    interpolator = precompile_NSE_3_5_7_soliton(_ε₂, _ε₃, 0.0, 0.0, L; use_M_values=use_M_values, M₀=M0, M₁=M1, warn_ignore=warn_ignore)
    _possible_solution_3_5_7 = (x) -> NSE_3_5_7_soliton(x, 0.0, 0.0, 0.0, 0.0, 0.0, interpolator)

    # x-coordinate correctrion
    shifted_grid_solution_3_5_7 = abs.(_possible_solution_3_5_7.(x_grid))
    bottom_gap_shifted = find_threshold(shifted_grid_solution_3_5_7, N_interpolation_points)
    x_reduced, y_reduced, _ = shift_pulse_to_center(
        x_grid,
        shifted_grid_solution_3_5_7;
        y_threshold=maximum(shifted_grid_solution_3_5_7) - bottom_gap_shifted
    )
    x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
    spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
    minimization_function = xi -> -spline_interpolator(xi)
    result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
    x_max_shifted, _ = result.minimizer, -result.minimum
    _z_0 = -(x_max_shifted - x_max)

    debug_flag && println("Солитон построен, но теперь его интерполированный максимум \
        находится в точке x=$x_max_shifted. Переносим на $_z_0")

    interpolator = precompile_NSE_3_5_7_soliton(_ε₂, _ε₃, 0.0, 0.0, L; use_M_values=use_M_values, M₀=M0, M₁=M1, warn_ignore=warn_ignore)
    possible_solution_3_5_7 = (x) -> NSE_3_5_7_soliton(x, 0.0, 0.0, 0.0, 0.0, _z_0, interpolator)
    deshifted_grid_solution_3_5_7 = abs.(possible_solution_3_5_7.(x_grid))

    if debug_flag
        # Проверяем что передвинули куда надо
        bottom_gap_deshifted = find_threshold(deshifted_grid_solution_3_5_7, N_interpolation_points)
        x_reduced, y_reduced, _ = shift_pulse_to_center(
            x_grid,
            deshifted_grid_solution_3_5_7;
            y_threshold=maximum(deshifted_grid_solution_3_5_7) - bottom_gap_deshifted
        )
        x_range = round(x_reduced[1], digits=10):h:round(x_reduced[end], digits=10)
        spline_interpolator = cubic_spline_interpolation((x_range), y_reduced)
        minimization_function = xi -> -spline_interpolator(xi)
        result = optimize(minimization_function, minimum(x_reduced), maximum(x_reduced))
        x_max_deshifted, _ = result.minimizer, -result.minimum
        println("Солитон перестроен, и теперь его интерполированный максимум находится в точке x=$x_max_deshifted")
        println("Для верности, индекс максимума: $(argmax(deshifted_grid_solution_3_5_7))")
    end

    return circshift(possible_solution_3_5_7.(x_grid), argmax(abs_U) - argmax(deshifted_grid_solution_3_5_7))
    #-argmax(deshifted_grid_solution_3_5_7)
end
function argmin_observed_change(old_index, vector; debug_flag=false)
    new_index = argmin(vector)
    debug_flag && println("был индекс $old_index, стал индекс $new_index")
    return new_index, old_index == new_index
end
function construct_best_approximate_NSE_3_5_7_solution(
    x_grid::Vector{<:Real},
    h::Real,
    abs_U::Vector{<:Real},
    L::Real,
    M0_initial;
    M0_gap=1.0,
    M0_h=0.5,
    debug_flag=false,
)
    y_max = maximum(abs_U)
    M0_vector_length = div(M0_gap, M0_h)
    construction_function = (M0, M1) -> abs.(construct_approximate_NSE_3_5_7_solution(
        x_grid, h, abs_U, 0.0, L;
        use_M_values=true, M0=M0, use_M1_value=true, M1=M1, warn_ignore=true,
    ))
    error_function = soliton -> relative_error_to_amplitude(abs_U, soliton)
    middle_index = Int(div(M0_gap, M0_h) + 1)

    error_argmin = 0
    solutions_with_same_amplitude = Vector{Float64}([])

    for _ in 1:6 # Целевое количество шагов дихотомии
        stop_flag = false
        while ~stop_flag
            M0_vector = round.(collect(M0_initial-M0_vector_length*M0_h:M0_h:M0_initial+M0_vector_length*M0_h), digits=10)

            debug_flag && println("Массив поиска: ", M0_vector)

            M0_vector = filter(x -> x < 0, M0_vector)
            M1_vector = 4 .* (y_max^2 .- M0_vector)
            solutions_with_same_amplitude = construction_function.(M0_vector, M1_vector)
            errors = error_function.(solutions_with_same_amplitude)
            error_argmin, stop_flag = argmin_observed_change(middle_index, errors; debug_flag)

            debug_flag && println(errors)
            stop_flag || (
                M0_initial = M0_vector[error_argmin];
                debug_flag && println("Останавливатся рано. Перемещаемся на $M0_initial")
            )
            stop_flag && (
                M0_initial = M0_vector[error_argmin];
                M0_h /= M0_vector_length;
                debug_flag && println("Достигли равновесия, уменьшаем шаг. Теперь M0_initial=$M0_initial, M0_h=$M0_h")
            )
        end
    end
    return solutions_with_same_amplitude[error_argmin]
end
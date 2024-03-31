function construct_approximate_NSE_5_solution(
    x,
    U,
    ε_2,
    theta_0,
    z_0,
)
    y_max=maximum(abs.(U))
    possible_μ = 2/3 * (2 * ε_2 * y_max^4 + 3 * y_max^2)
    possible_solution_5 = NSE_5_soliton.(x, 0, 0, possible_μ/4, ε_2, theta_0, z_0)

    # x-coordinate correctrion
    _z_0 = x[argmax(abs.(possible_solution_5))]-x[argmax(abs.(U))]
    return NSE_5_soliton.(x, 0, 0, possible_μ/4, ε_2, theta_0, -_z_0)
end
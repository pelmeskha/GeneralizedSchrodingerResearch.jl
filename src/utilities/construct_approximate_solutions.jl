function construct_approximate_NSE_3_5_solution(
    x,
    t,
    U::Union{Vector{ComplexF64},Vector{<:Real}},
    k::Real,
    ε₂::Real,
    theta_0::Real,
    z_0::Real,
    L::Real,
)
    y_max = maximum(abs.(U))
    possible_μ = 2 / 3 * (2 * ε₂ * y_max^4 + 3 * y_max^2)
    _possible_ω = possible_μ / 4
    possible_ω = _possible_ω + k^2
    possible_solution_5 = NSE_3_5_soliton.(x, t, k, possible_ω, ε₂, theta_0, z_0, cycle=true, L=L)

    # x-coordinate correctrion
    _z_0 = x[argmax(abs.(possible_solution_5))] - x[argmax(abs.(U))]
    return NSE_3_5_soliton.(x, t, k, possible_ω, ε₂, theta_0, -_z_0, cycle=true, L=L)
end
"""
  #TODO:
  - Сделать поддержку t с зацикливанием;
  - Проверить формулы для M0, M1 - нужный максимум достигается не так точно, как хотелось бы
    возможно стоит более просто подбирать эпсилоны по соотношению для максимума.
"""
function construct_approximate_NSE_3_5_7_solution(
    x,
    t,
    U::Union{Vector{ComplexF64},Vector{<:Real}},
    k::Real,
    ε₂::Real,
    theta_0::Real,
    z_0::Real,
    ξ₀::Real,
    L::Real,
)
    y_max = maximum(abs.(U))

    M0 = -(sqrt(4.0 * ε₂^2 * y_max^2 + 6.0 * ε₂ * y_max + 9.0) + 2 * ε₂ + 3.0) / (2.0 * ε₂)
    M1 = (2.0 * (3.0 + 4.0 * ε₂ * y_max + sqrt(9.0 + 6.0 * ε₂ * y_max + 4.0 * ε₂^2 * y_max^2))) / ε₂
    _ε₂, _ε₃ = M0_M1_to_ε2_ε3(M0, M1)
    println("Параметры подобраны для ε₂ = ", _ε₂, ", ε₃ = ", _ε₃)

    #ω = - (k^2 - 1.0 / 12.0 * (M0 * M1) / (M1 + 6 * M0) - 1.0 / 6.0 * M0)
    interpolator = precompile_NSE_3_5_7_soliton(_ε₂, _ε₃, z_0, ξ₀, L)
    _possible_solution_5_7 = (x) -> NSE_3_5_7_soliton(x, 0.0, 0.0, 0.0, theta_0, z_0, interpolator)
    possible_solution_5_7 = _possible_solution_5_7.(x)
    # x-coordinate correctrion
    _z_0 = x[argmax(abs.(possible_solution_5_7))] - x[argmax(abs.(U))]
    println("z0 = ",_z_0)
    interpolator = precompile_NSE_3_5_7_soliton(_ε₂, _ε₃, 0.0, ξ₀, L)
    _possible_solution_5_7 = (x) -> NSE_3_5_7_soliton(x, 0.0, 0.0, 0.0, theta_0, -_z_0, interpolator)
    return _possible_solution_5_7.(x)
end
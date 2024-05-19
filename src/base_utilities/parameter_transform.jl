function M0_M1_to_ε2_ε3(M0::Real, M1::Real)
    ε2 = 3.0 / (4.0 * M0) * (M1 / (M1 + 6 * M0) - 2.0)
    ε3 = 4.0 / (M0 * (M1 + 6 * M0))
    return (ε2, ε3)
end
function ε2_ε3_to_M0_M1(ε2::Real, ε3::Real)
    M0_a = (-4*ε2-sqrt(2)*sqrt(8*ε2^2-27*ε3))/(9*ε3)
    M1_a = (4*sqrt(2)*sqrt(8*ε2^2-27*ε3))/(3*ε3)

    M0_b = (-4*ε2+sqrt(2)*sqrt(8*ε2^2-27*ε3))/(9*ε3)
    M1_b = -(4*sqrt(2)*sqrt(8*ε2^2-27*ε3))/(3*ε3)

    (M0, M1) = M1_a > 0.0 ? (M0_a, M1_a) : (M0_b, M1_b)
    return (M0, M1)
end
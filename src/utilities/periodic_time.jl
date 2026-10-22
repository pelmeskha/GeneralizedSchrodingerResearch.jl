"""
    periodic_time(z, T, Tspan, c0, w0)

    Map T to the periodic image closest to the analytical pulse center.

    The period is the FFT domain length:
        L_T = Tspan[2] - Tspan[1]

    Requires c0 != 0.
"""
function periodic_time(z, T, Tspan, c0, w0)
    L_T = Tspan[2] - Tspan[1]

    L_T > 0 || throw(ArgumentError("Tspan must have positive length"))
    c0 != 0 || throw(ArgumentError("c0 must be nonzero"))

    center = (z - w0) / c0

    return center + mod(T - center + L_T/2, L_T) - L_T/2
end

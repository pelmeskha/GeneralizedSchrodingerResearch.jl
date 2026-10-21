using GeneralizedSchrodingerResearch.Solvers: solve_coupled_gl
include("plotting.jl")
using .CoupledGLPlotting
using GeneralizedSchrodingerResearch.AnalyticalSolutions.CoupledGLAnalytical

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
"""
    coupled_gl_residual(
        u_solution, v_solution, z, t;
        Δβ, δ, β2, g, Ω_g, γ,
        hz=1e-3, ht=1e-3
    )

    Evaluate the coupled Ginzburg–Landau PDE residuals at (z, t).

    Solutions must be callable as:
        u_solution(z, t)
        v_solution(z, t)

    Derivatives use fourth-order centered finite differences.
    The functions are evaluated at z ± hz, z ± 2hz,
    and t ± ht, t ± 2ht.

    Returns:
        residual_u, residual_v : complex PDE residuals
        error_u, error_v       : absolute residual magnitudes

    A smooth exact solution should produce residuals near zero,
    up to differentiation and floating-point errors.
"""
function coupled_gl_residual(
    u_solution, v_solution;
    z, t, Δβ, δ, β2, g, Ω_g, γ,
    hz=1e-3,
    ht=1e-3
)
    hz > 0 || throw(ArgumentError("hz must be positive"))
    ht > 0 || throw(ArgumentError("ht must be positive"))
    Ω_g != 0 || throw(ArgumentError("Ω_g must be nonzero"))

    function field_derivatives(field)
        f = field(z, t)

        fz_minus2 = field(z - 2*hz, t)
        fz_minus1 = field(z - hz, t)
        fz_plus1  = field(z + hz, t)
        fz_plus2  = field(z + 2*hz, t)

        ft_minus2 = field(z, t - 2*ht)
        ft_minus1 = field(z, t - ht)
        ft_plus1  = field(z, t + ht)
        ft_plus2  = field(z, t + 2*ht)

        f_z = (
            fz_minus2 - 8*fz_minus1 +
            8*fz_plus1 - fz_plus2
        ) / (12*hz)

        f_t = (
            ft_minus2 - 8*ft_minus1 +
            8*ft_plus1 - ft_plus2
        ) / (12*ht)

        f_tt = (
            -ft_plus2 + 16*ft_plus1 - 30*f +
            16*ft_minus1 - ft_minus2
        ) / (12*ht^2)

        return f, f_z, f_t, f_tt
    end

    u, u_z, u_t, u_tt = field_derivatives(u_solution)
    v, v_z, v_t, v_tt = field_derivatives(v_solution)

    residual_u = (
        u_z
        - 1im * Δβ/2 * u
        + δ * u_t
        + 1im * β2/2 * u_tt
        - g/2 * u
        - g/(2*Ω_g^2) * u_tt
        - 1im * γ * (abs2(u) + (2/3)*abs2(v)) * u
        - 1im * γ/3 * v^2 * conj(u)
    )

    residual_v = (
        v_z
        - 1im * Δβ/2 * v
        + δ * v_t
        + 1im * β2/2 * v_tt
        - g/2 * v
        - g/(2*Ω_g^2) * v_tt
        - 1im * γ * (abs2(v) + (2/3)*abs2(u)) * v
        - 1im * γ/3 * u^2 * conj(v)
    )

    return (
        residual_u=residual_u,
        residual_v=residual_v,
        error_u=abs(residual_u),
        error_v=abs(residual_v)
    )
end
"""
    Evaluate residual magnitudes on a grid.

    Output matrices have shape:
        (length(t_values), length(z_values))
    # Check if equations satisfied
    parameters = (Δβ=Δβ, δ=δ, β2=β2, g=g, Ω_g=Ω_g, γ=γ)

    check = coupled_gl_residual_grid(
        true_u,
        true_v;
        z_values = range(0.0, 40.0, length=101),
        t_values = range(Tspan[1], Tspan[2], length=256),
        parameters...,
        hz=1e-4,
        ht=1e-4
    )

    println("True solution check:")
    println("Maximum |R_u|: ", check.max_error_u)
    println("Maximum |R_v|: ", check.max_error_v)
"""
function coupled_gl_residual_grid(
    u_solution, v_solution;
    z_values, t_values,
    kwargs...
)
    error_u = Matrix{Float64}(
        undef, length(t_values), length(z_values)
    )
    error_v = similar(error_u)

    for (i, zi) in enumerate(z_values)
        for (j, tj) in enumerate(t_values)
            result = coupled_gl_residual(
                u_solution, v_solution;
                z=zi, t=tj,
                kwargs...
            )

            error_u[j, i] = result.error_u
            error_v[j, i] = result.error_v
        end
    end

    return (
        error_u=error_u,
        error_v=error_v,
        max_error_u=maximum(error_u),
        max_error_v=maximum(error_v)
    )
end

# domain
Tspan = (-100.0, 100.0)
L = Tspan[2]-Tspan[1]
zspan = (0.0, 300.0)
dz    = 0.025
N_T   = 256

# pulse parameters
k = 12.6
c0 = 1.4
mu = 0.97
I1 = 0.1
ω = 0
w0 = 0 # displacement

# equation parameters. Construct corresponding equation
β2 = 0.025; # should be ~ 0.025
Ω_g = 10.0; # should be ~ [10, 100]

# from restrictions
δ  = (1 - β2*c0*ω) / c0;                                  # ~ [0, 0.1]. More -> faster | β2 c0    ω
g  = 0;                                                   # ~ [0.01, 0.5]
γ  = mu^4 * c0^2 * β2 / (32*I1)                           # ~ 0.0013                   | β2 c0 mu    I1
Δβ = (-β2*c0^3*mu^2 + 2*β2*c0*ω^2 + 4*c0*k - 4*ω)/(2*c0); # ~ [0, 10]                  | β2 c0 mu ω     k

# helper
# A = 2 * sqrt(2 * I1) / abs(mu)
xi = coupled_gl_xi(; Δβ, β2, γ, c0, k, ω)
sn_pars = sn_parameters(xi.xi1, xi.xi2, c0 * L)

u0(T) = sn_solution(0.0, T, sn_pars; c0, w0, k, ω)
# u0(T) = A * tanh(mu/2 * (0 - c0*T - w0)) * exp(1im * (- ω*T)) - A * tanh(mu/2 * (0 - c0*T - w0 - c0 * L/2)) * exp(1im * (- ω*T))
v0(T) = 0 # sin(T) * exp(-(-c0*T-w0)^2) # 

true_u(z, T) = sn_solution(z, T, sn_pars; c0, w0, k, ω)
# true_u(z, T) = A * tanh(mu * (z - c0*T - w0) / 2) * exp(1im * (k*z - ω*T))
# true_u(z, T) = A_per * jacobi_sn(kappa_per * (z - c0*T - w0), kp_per) * exp(1im * (k*z - ω*T))
true_v(z, T) = 0
# Periodic analytical reference.
true_u_periodic(z, T) = true_u(z, periodic_time(z, T, Tspan, c0, w0))
true_v_periodic(z, T) = true_v(z, periodic_time(z, T, Tspan, c0, w0))

# δ  = 20
# γ  = 0.0013
# Δβ = 5
println("β2 = $β2, Ω_g = $Ω_g, Δβ = $Δβ, γ = $γ, δ = $δ")


println("Примерный критерий сходимости: $(dz/( (Tspan[2]-Tspan[1]) /N_T)^2) << $(0.2/(abs(β2)*pi^2)) ?")

T, z, u, v = solve_coupled_gl(zspan, Tspan, dz, N_T, u0, v0, Δβ, δ, β2, g, Ω_g, γ)

u_metrics = field_error(u, T, z, true_u_periodic)
println(
    "u vs analytical: max modulus error = ", round(u_metrics.max_modulus_error, sigdigits=3),
    ", max complex error = ", round(u_metrics.max_complex_error, sigdigits=3),
    ", max relative error = ", round(u_metrics.max_relative_error_percent, sigdigits=3), " %"
)

plot_coupled_gl(T, z, u, v; prefix="2D plot actual", cmap=:plasma)
plot_coupled_gl(T, z, true_u_periodic, true_v_periodic; prefix="2D plot analytical", cmap=:plasma) # need other function?
plot_profiles(
    T_grid = T,
    z_grid = z,
    z_timestamps = [zspan[1], zspan[2]],
    u_solutions = [u, true_u_periodic],
    # v_solutions = [v, true_v_periodic], # commented = empty array by default. No plots for v component
    prefixes=["actual", "analytical"],
)
plot_true_errors(T, z, u, v, true_u_periodic, true_v_periodic; prefix="error plots")
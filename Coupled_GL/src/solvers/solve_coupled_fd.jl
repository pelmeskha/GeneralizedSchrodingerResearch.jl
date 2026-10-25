using LinearAlgebra: Tridiagonal, lu, ldiv!, mul!
using ProgressMeter

# Finite-difference solver in T with non-periodic boundaries.
#
# Equation: ∂z u = L u + N(u, v),  ∂z v = L v + N(v, u) with
#   L = g/2 + iΔβ/2 - δ ∂T + (g/(2Ω_g²) - iβ2/2) ∂T²,
# the same operator as D̃(ω) in solve_coupled_gl with ω ↔ -i∂T.
#
# Spatial scheme: second-order central differences on the inclusive grid T_1..T_N.
# Boundary conditions:
#   :dirichlet  (1st kind) u = v = 0 at both ends (pinned nodes).
#   :neumann    (2nd kind) ∂T u = ∂T v = 0 at both ends (mirror ghost nodes).
# Splitting: Strang, L(Δz/2) N(Δz) L(Δz/2). The linear half-steps use Crank–Nicolson
# with step Δz/2, the nonlinear step uses RK4 as in solve_coupled_gl.
function solve_coupled_gl_fd(
    zspan::Tuple{Real,Real},     # (zmin, zmax)
    Tspan::Tuple{Real,Real},     # (Tmin, Tmax), both ends are grid nodes
    dz::Real,                    # step size in z
    N_T::Int,                    # number of T nodes including both ends (≥ 3)
    u0::Function,                # u(0,T)
    v0::Function,                # v(0,T)
    Δβ::Real, δ::Real, β2::Real, # dispersion / birefringence
    g::Real, Ω_g::Real,          # gain parameters
    γ::Real;                     # nonlinear coefficient
    bc::Symbol = :dirichlet,     # :dirichlet or :neumann
    show_progress = true,
    store = true,                # false: do not keep the fields
    store_every::Int = 1,        # keep every store_every-th z plane (returned z is subsampled to match)
    on_step = nothing            # called as on_step(k, z[k], u, v) after every z step
)
    bc in (:dirichlet, :neumann) || throw(ArgumentError("bc must be :dirichlet or :neumann, got :$bc"))
    store_every >= 1 || throw(ArgumentError("store_every must be at least 1, got $store_every"))
    N_T >= 3 || throw(ArgumentError("N_T must be at least 3, got $N_T"))

    # ---- T grid (inclusive) ----
    T  = collect(range(Tspan[1], Tspan[2], length=N_T))
    h  = (Tspan[2] - Tspan[1]) / (N_T - 1)

    # ---- z grid ----
    N_z = round(Int, (zspan[2] - zspan[1]) / dz) + 1
    z   = range(zspan[1], zspan[2], length=N_z)
    Δz  = step(z)
    stored = 1:store_every:N_z
    z_out = z[stored]

    # ---- linear operator, tridiagonal in T ----
    lower, diag, upper = fd_linear_bands(N_T, h, Δβ, δ, β2, g, Ω_g, bc)

    # ---- Crank–Nicolson half-steps: (I - τL) u⁺ = (I + τL) u, τ = Δz/4 for the step Δz/2 ----
    τ = Δz / 4
    A_minus = Tridiagonal(-τ .* lower, 1 .- τ .* diag, -τ .* upper)
    A_plus  = Tridiagonal( τ .* lower, 1 .+ τ .* diag,  τ .* upper)
    F_lin   = lu(A_minus)
    buf     = Vector{ComplexF64}(undef, N_T)

    # ---- initial condition ----
    u = ComplexF64.(u0.(T))
    v = ComplexF64.(v0.(T))
    if bc === :dirichlet
        u[1] = u[end] = 0
        v[1] = v[end] = 0
    end

    # ---- output arrays ----
    U_out = store ? zeros(ComplexF64, N_T, length(stored)) : nothing
    V_out = store ? zeros(ComplexF64, N_T, length(stored)) : nothing
    if store
        U_out[:,1] .= u
        V_out[:,1] .= v
    end
    isnothing(on_step) || on_step(1, z[1], u, v)

    # ---- buffers for RK4 ----
    f1u = similar(u); f2u = similar(u); f3u = similar(u); f4u = similar(u)
    f1v = similar(v); f2v = similar(v); f3v = similar(v); f4v = similar(v)
    utmp = similar(u); vtmp = similar(v)

    function nonlinear_rhs!(du, dv, u_cur, v_cur)
        @. du = 1im*γ * (
            (abs2(u_cur) + (2/3)*abs2(v_cur)) * u_cur +
            (1/3) * v_cur^2 * conj(u_cur)
        )
        @. dv = 1im*γ * (
            (abs2(v_cur) + (2/3)*abs2(u_cur)) * v_cur +
            (1/3) * u_cur^2 * conj(v_cur)
        )
    end

    function rk4_nonlinear!(u, v, dz_local)
        nonlinear_rhs!(f1u, f1v, u, v)
        @. utmp = u + 0.5*dz_local*f1u
        @. vtmp = v + 0.5*dz_local*f1v
        nonlinear_rhs!(f2u, f2v, utmp, vtmp)
        @. utmp = u + 0.5*dz_local*f2u
        @. vtmp = v + 0.5*dz_local*f2v
        nonlinear_rhs!(f3u, f3v, utmp, vtmp)
        @. utmp = u + dz_local*f3u
        @. vtmp = v + dz_local*f3v
        nonlinear_rhs!(f4u, f4v, utmp, vtmp)
        @. u += (dz_local/6) * (f1u + 2*f2u + 2*f3u + f4u)
        @. v += (dz_local/6) * (f1v + 2*f2v + 2*f3v + f4v)
    end

    function linear_half_step!(w)
        mul!(buf, A_plus, w)
        ldiv!(w, F_lin, buf)
    end

    # ---- main loop with symmetric splitting ----
    iter = 2:N_z
    p = show_progress ? Progress(length(iter); desc="Propagating (FD)") : nothing

    for k in iter
        linear_half_step!(u)
        linear_half_step!(v)

        rk4_nonlinear!(u, v, Δz)

        linear_half_step!(u)
        linear_half_step!(v)

        if store && (k - 1) % store_every == 0
            idx = (k - 1) ÷ store_every + 1
            U_out[:,idx] .= u
            V_out[:,idx] .= v
        end
        isnothing(on_step) || on_step(k, z[k], u, v)

        show_progress && next!(p)
    end

    return T, z_out, U_out, V_out
end

# Bands of L (sub-, main, super-diagonal) for the central-difference operator on N nodes.
function fd_linear_bands(N, h, Δβ, δ, β2, g, Ω_g, bc)
    c  = g/2 + 1im*Δβ/2
    a1 = -δ
    a2 = g/(2*Ω_g^2) - 1im*β2/2

    lower = fill(-a1/(2h) + a2/h^2, N-1)
    diag  = fill(c - 2*a2/h^2, N)
    upper = fill( a1/(2h) + a2/h^2, N-1)

    if bc === :dirichlet
        diag[1] = diag[end] = 0
        upper[1] = 0
        lower[end] = 0
    else
        upper[1] = 2*a2/h^2
        lower[end] = 2*a2/h^2
    end
    return ComplexF64.(lower), ComplexF64.(diag), ComplexF64.(upper)
end

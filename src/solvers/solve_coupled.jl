using FFTW
using ProgressMeter

function solve_coupled_gl(
    zspan::Tuple{Real,Real},      # (zmin, zmax)
    Tspan::Tuple{Real,Real},      # (Tmin, Tmax)
    dz::Real,                     # step size in z
    N_T::Int,                     # number of T‑grid points (even)
    u0::Function,                 # u(0,T)
    v0::Function,                 # v(0,T)
    Δβ::Real, δ::Real, β2::Real,  # dispersion / birefringence
    g::Real, Ω_g::Real,           # gain parameters
    γ::Real;                      # nonlinear coefficient
    show_progress = true
)
    # ---- T grid and frequency grid ----
    L_T = Tspan[2] - Tspan[1]
    dT  = L_T / N_T
    T   = collect(range(Tspan[1], step=dT, length=N_T))   # periodic grid
    ω   = 2π .* fftfreq(N_T, dT)                          # FFTW frequency ordering

    # ---- z grid ----
    N_z = round(Int, (zspan[2] - zspan[1]) / dz) + 1
    z   = range(zspan[1], zspan[2], length=N_z)
    Δz  = step(z)                                          # actual step used

    # ---- linear operator in frequency domain ----
    D̃(ω_) = g/2 - g*ω_^2/(2*Ω_g^2) + 1im*(Δβ/2 - δ*ω_ + β2/2*ω_^2)
    Lin_half = exp.(D̃.(ω) .* Δz/2)

    # ---- initial condition ----
    u = ComplexF64.(u0.(T))
    v = ComplexF64.(v0.(T))

    # ---- output arrays ----
    U_out = zeros(ComplexF64, N_T, N_z)
    V_out = zeros(ComplexF64, N_T, N_z)
    U_out[:,1] .= u
    V_out[:,1] .= v

    # ---- FFT plans (create proper plans for in-place and out-of-place transforms) ----
    # For in-place transformation of u and v separately
    P_u    = plan_fft!(u)        # in-place FFT for u
    Pinv_u = plan_ifft!(u)       # in-place inverse FFT for u
    P_v    = plan_fft!(v)        # in-place FFT for v  
    Pinv_v = plan_ifft!(v)       # in-place inverse FFT for v

    # ---- buffers for RK4 ----
    f1u = similar(u); f2u = similar(u); f3u = similar(u); f4u = similar(u)
    f1v = similar(v); f2v = similar(v); f3v = similar(v); f4v = similar(v)
    utmp = similar(u); vtmp = similar(v)

    # nonlinear right‑hand side (in‑place, writes into du, dv)
    function nonlinear_rhs!(du, dv, u_cur, v_cur)
        @. du = 1im*γ * ( (abs2(u_cur) + (2/3)*abs2(v_cur)) * u_cur +
                          (1/3) * v_cur^2 * conj(u_cur) )
        @. dv = 1im*γ * ( (abs2(v_cur) + (2/3)*abs2(u_cur)) * v_cur +
                          (1/3) * u_cur^2 * conj(v_cur) )
    end

    # RK4 step for the nonlinear ODE (in‑place update of u, v)
    function rk4_nonlinear!(u, v, dz_local)
        # k1
        nonlinear_rhs!(f1u, f1v, u, v)
        # k2
        @. utmp = u + 0.5*dz_local*f1u
        @. vtmp = v + 0.5*dz_local*f1v
        nonlinear_rhs!(f2u, f2v, utmp, vtmp)
        # k3
        @. utmp = u + 0.5*dz_local*f2u
        @. vtmp = v + 0.5*dz_local*f2v
        nonlinear_rhs!(f3u, f3v, utmp, vtmp)
        # k4
        @. utmp = u + dz_local*f3u
        @. vtmp = v + dz_local*f3v
        nonlinear_rhs!(f4u, f4v, utmp, vtmp)
        # assemble
        @. u += (dz_local/6) * (f1u + 2*f2u + 2*f3u + f4u)
        @. v += (dz_local/6) * (f1v + 2*f2v + 2*f3v + f4v)
    end

    # ---- main loop with symmetric splitting ----
    iter = 2:N_z
    p = show_progress ? Progress(length(iter); desc="Propagating") : nothing

    for k in iter
        # --- half linear step ---
        # Apply in-place: FFT, multiply, inverse FFT
        mul!(u, P_u, u)          # u = fft(u)
        u .*= Lin_half
        mul!(u, Pinv_u, u)       # u = ifft(u)
        
        mul!(v, P_v, v)          # v = fft(v)
        v .*= Lin_half
        mul!(v, Pinv_v, v)       # v = ifft(v)

        # --- full nonlinear step ---
        rk4_nonlinear!(u, v, Δz)

        # --- second half linear step ---
        mul!(u, P_u, u)          # u = fft(u)
        u .*= Lin_half
        mul!(u, Pinv_u, u)       # u = ifft(u)
        
        mul!(v, P_v, v)          # v = fft(v)
        v .*= Lin_half
        mul!(v, Pinv_v, v)       # v = ifft(v)

        # store
        U_out[:,k] .= u
        V_out[:,k] .= v

        show_progress && next!(p)
    end

    return T, z, U_out, V_out
end
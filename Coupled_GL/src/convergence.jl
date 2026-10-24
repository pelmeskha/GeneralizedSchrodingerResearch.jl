function max_relative_error_over_run(setup, dz, N_T)
    (; Tspan, zspan, args, u0, v0, true_u_periodic) = setup
    dT = (Tspan[2] - Tspan[1]) / N_T
    T = collect(range(Tspan[1], step=dT, length=N_T))
    max_err = 0.0
    max_ref = 0.0
    check_every = max(1, round(Int, 0.01 / dz))

    function check_step!(k, z_k, u, v)
        (k - 1) % check_every == 0 || return
        for j in eachindex(T)
            ref = true_u_periodic(z_k, T[j])
            max_err = max(max_err, abs(u[j] - ref))
            max_ref = max(max_ref, abs(ref))
        end
    end

    t0 = time()
    _, z, _, _ = solve_coupled_gl(
        zspan, Tspan, dz, N_T, u0, v0,
        args.Δβ, args.δ, args.β2, args.g, args.Ω_g, args.γ;
        show_progress=false, store=false, on_step=check_step!,
    )
    runtime = time() - t0

    return (
        dT=dT,
        N_z=length(z),
        max_relative_error_percent=100 * max_err / max_ref,
        runtime=runtime,
    )
end

function save_convergence_csv(path, header, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, header)
        for r in rows
            println(io, join(r, ","))
        end
    end
end

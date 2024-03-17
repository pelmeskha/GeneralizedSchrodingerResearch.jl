#= 
    По начальному условию и заданным промежуткам решает начально-краевую задачу методом Фурье
=#
function fourier_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function;
    ε_2 = 1,
    ε_3 = 1,
    # filtration parameters
    filtration_flag::Bool = false,
    filtration_step::Int = 10,
    filtration_factor::Float64 = 0.5,
    l_nominal=100,
    # tolerance calculations
    tolerance_flag = false,
    analytical_solution = [],
)
    theta = 0.5
    L = xspan[2] - xspan[1]
    T = tspan[2] - tspan[1]
    N = Int(L / h)
    N_x = N + 1
    N_t = Int(round(T / tau + 1, digits = 1))

    j = range(-N / 2, stop = N / 2 , length = N_x)
    x = collect(j .* h)[1:end-1]
    t = range(tspan[1], tspan[2], length = N_t)

    mun = collect(2 * pi / L .* range(-N / 2, stop = N / 2 - 1, length = N))
    direct = (h / L .* exp.(-1im .* mun * x'))
    inverse = exp.(1im .* mun * x')
    M=exp.(-1im.*mun.^2 .*tau)

    U = initial_function.(x)

    (sum(abs.(U) .> abs(U[1])) * sum(abs.(U) .> abs(U[end])) != 0) ||
        throw(AssertionError("Initial pulse out of x interval."))
    
    ε_tail = 1e-5
    (abs(U[1]) < ε_tail) & (abs(U[end]) < ε_tail)  ||
        throw(AssertionError("Pulse tail modulus exceeds $ε_tail with the value of $(max(abs(U[1]), abs(U[end]))). Consider larger x-interval."))

    if filtration_flag
        power_dissipated=0.0
    end
    if tolerance_flag
        tolerance = zeros(N_t)
        tolerance[1] = maximum(abs.(U)-abs.(analytical_solution.(x,tspan[1])))
    end
    @showprogress for i in 1:N_t-1
        if filtration_flag
            if mod(i*tau, filtration_step) == 0.0
                U, power = filtration(
                    U,
                    x,
                    filtration_factor,
                    l_nominal,
                )
                power_dissipated += power
            end
        end

        U=inverse*(
            M.*(
                direct*(
                    exp.(
                        1im*tau* ((abs.(U)).^2 + ε_2*(abs.(U)).^4 + ε_3*(abs.(U)).^6 ) 
                    ).*U
                )
            )
        )
        if tolerance_flag
            # Percentage tolerance
            tolerance[i+1] = maximum( 
                (abs.(analytical_solution.(x,i*tau)) - abs.(U))
            ) / maximum(abs.(analytical_solution.(x,i*tau))) * 100
        end
    end
    return(
        x,
        U,
        filtration_flag ? power_dissipated : nothing,
        tolerance_flag ? tolerance : nothing,
    )
end
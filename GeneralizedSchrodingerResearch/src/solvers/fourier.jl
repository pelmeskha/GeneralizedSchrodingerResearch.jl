#= 
    По начальному условию и заданным промежуткам решает начально-краевую задачу методом Фурье
=#
function fourier_solve(
    tspan::Tuple{Float64, Float64},
    xspan::Tuple{Float64, Float64},
    tau::Float64,
    h::Float64,
    initial_function;
    ε_2::Float64 = 0.0,
    ε_3::Float64 = 0.0,
    # filtration parameters
    filtration_flag::Bool = false,
    filtration_time::Float64 = 10.0,
    filtration_factor::Float64 = 0.5,
    l_nominal::Float64=100.0,
    # tolerance calculations
    tolerance_flag = false,
    analytical_solution = [],
    # record integrals
    integrals_flag = false,
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
        t_slider=filtration_time
        I1_dissipated=zeros(N_t)
        I2_dissipated=zeros(N_t)
    end
    if tolerance_flag
        tolerance = zeros(N_t)
    end
    if integrals_flag
        I_1 = zeros(N_t)
        I_2 = zeros(N_t)
        I_1[1]=integral_1(U,h)
        I_2[1]=integral_2(U,h)
    end
    @showprogress for i in 1:N_t
        if tolerance_flag
            if isa(analytical_solution,Real) # calculate the "pike" tolerance
                tolerance[i]=(maximum(abs.(U))-analytical_solution)/maximum(abs.(U)) * 100
            else
                # Percentage tolerance
                tolerance[i] = maximum( 
                    (abs.(analytical_solution.(x,i*tau)) - abs.(U))
                ) / maximum(abs.(analytical_solution.(x,i*tau))) * 100
            end
        end
        if integrals_flag
            I_1[i] = integral_1(U,h)
            I_2[i] = integral_2(U,h)
        end
        if filtration_flag
            if i*tau >= t_slider
                t_slider+=filtration_time
                U, (power_I1, power_I2) = filtration(
                    U,
                    x,
                    filtration_factor,
                    l_nominal,
                )
                I1_dissipated[i] = power_I1
                I2_dissipated[i] = power_I2
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
    end
    return(
        x,
        t,
        U,
        filtration_flag ? (cumsum(I1_dissipated), cumsum(I2_dissipated)) : (nothing, nothing),
        tolerance_flag ? tolerance : nothing,
        integrals_flag ? (I_1, I_2) : (nothing, nothing),
    )
end
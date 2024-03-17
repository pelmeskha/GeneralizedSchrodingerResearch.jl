#= 
    По начальному условию и заданным промежуткам решает начально-краевую задачу методом Фурье
=#
function fourier_solve(tspan,xspan,tau,h,initial_function; ε_2 = 1, ε_3 = 1)
    theta = 0.5
    L = xspan[2] - xspan[1]
    T = tspan[2] - tspan[1]
    N = Int(L / h)
    N_x = N + 1
    N_t = Int(T / tau + 1)

    j = range(-N / 2, stop = N / 2 , length = N_x)
    x = collect(j .* h)
    t = range(tspan[1], tspan[2], length = N_t)

    mun = collect(2 * pi / L .* range(-N / 2, stop = N / 2 - 1, length = N))
    direct = (h / L .* exp.(-1im .* mun * x[1:end-1]'))
    inverse = exp.(1im .* mun * x[1:end-1]')
    M=exp.(-1im.*mun.^2 .*tau)

    U = initial_function(x[1:end-1])
    @showprogress for i in 1:N_t-1
        if mod(i*tau, 2) == 0
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
    return(x, U)
end
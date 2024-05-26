"""
    По начальному условию и заданным промежуткам решает начально-краевую задачу Фурье или
    конечно-разностным методом. Реализация на CUDA.

    Опционально:
        - сравнивает численное решение с заданным аналитическим;
        - вычисляет интегралы в процессе моделирования;
        - реализует алгоритм фильтрации излучения;
        - сохраняет решение в заданные моменты времени.
    
    #TODO:
    - Совместить с основным solve.
"""
function cuda_solve(
    tspan::Tuple{Real, Real},
    xspan::Tuple{Real, Real},
    tau::Real,
    h::Real,
    initial_function;
    method::String = "fourier",
    ε_2::Real = 0.0,
    ε_3::Real = 0.0,
    # filtration parameters
    filtration_flag::Bool = false,
    filtration_time::Real = 10.0,
    filtration_factor::Union{Real, Function} = 1.0,
    filtration_end_t::Real = tspan[2],
    l_nominal::Real=100.0,
    # pulse_maximum calculations
    pulse_maximum_flag = false,
    # record integrals
    integrals_flag = false,
    # times of interest
    capture_times = [],
    # live plot
    live_plot_solution_flag = false,
    live_plot_tau = 1.0,
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

    if method == "fourier"
        mun = collect(2 * pi / L .* range(-N / 2, stop = N / 2 - 1, length = N))
        direct = (h / L .* exp.(-1im .* mun * x'))
        inverse = exp.(1im .* mun * x')
        fourier_ratio=exp.(-1im.*mun.^2 .*tau)
        M = inverse * (fourier_ratio .* direct)
    else # method == "finite_difference"
        S = create_FD_matrix(N)
        I = Diagonal(ones(N))
        r = tau/h^2
        M = (I-1im*r*theta*S)^-1 * (I+1im*r*(1-theta)*S)
    end

    U = initial_function.(x)

    (sum(abs.(U) .> abs(U[1])) * sum(abs.(U) .> abs(U[end])) != 0) ||
        throw(AssertionError("Initial pulse out of x interval."))
    
    ε_tail = 1e-5
    (abs(U[1]) < ε_tail) & (abs(U[end]) < ε_tail)  ||
    @warn "Pulse tail modulus exceeds $ε_tail with the value of $(max(abs(U[1]), abs(U[end]))). Consider larger x-interval."

    if filtration_flag
        t_slider=filtration_time
        I1_dissipated=zeros(N_t)
        I2_dissipated=zeros(N_t)
    end
    if pulse_maximum_flag
        pulse_maximum = zeros(N_t)
    end
    if integrals_flag
        I_1 = zeros(N_t)
        I_2 = zeros(N_t)
        I_1[1]=integral_1(U,h)
        I_2[1]=integral_2(U,h)
    end
    capture_times_flag=false
    if ~isempty(capture_times)
        (any(capture_times.>tspan[2]) || any(capture_times.<tspan[1])) && @error "capture_times contains an element out of modeling time"
        capture_times_flag=true
        push!(capture_times,Inf)
        capture_times_slider=1
        t_capture_entry=capture_times[capture_times_slider]
        U_set=[]
    end
    cuda_U = CuArray(U)
    cuda_M = CuArray(M)
    if live_plot_solution_flag
        observable_U = Observable(abs.(U))
        observable_x = Observable(x)
        live_plot_t_slider = 0.0
        live_plot = Figure()
        live_axis = Axis(live_plot[1, 1], limits = (minimum(x), maximum(x), 0.0, 1.0))
        lines!(live_axis, observable_x, observable_U)
        display(live_plot)
    end
    #@showprogress for i in 1:N_t
    for i in 1:N_t
        t_current=(i-1)*tau
        if live_plot_solution_flag && (t_current ≥ live_plot_t_slider)
            observable_U[] = abs.(Array(cuda_U))
            live_plot_t_slider += live_plot_tau
        end
        if capture_times_flag && t_capture_entry≠Inf
            if t_current ≥ t_capture_entry
                host_U=Array(cuda_U)
                push!(U_set, host_U)
                capture_times_slider+=1
                t_capture_entry=capture_times[capture_times_slider]
            end
        end
        if pulse_maximum_flag
            pulse_maximum[i]=cuda_maximum(cuda_U)
        end
        if filtration_flag && (t_current <= filtration_end_t)
            if t_current ≥ t_slider
                t_slider+=filtration_time
                cuda_U, (power_I1, power_I2) = cuda_filtration(
                    cuda_U,
                    h,
                    isa(filtration_factor, Function) ? filtration_factor(t_current) : filtration_factor,
                    l_nominal,
                )
                I1_dissipated[i] = power_I1
                I2_dissipated[i] = power_I2
            end
        end

        cuda_V = cuda_calculate_V(cuda_U, tau, ε_2, ε_3)
        cuda_U = cuda_matrix_vector_multiplication(cuda_M,cuda_V)

        if integrals_flag
            I_1[i] = cuda_integral_1(cuda_U,h)
            I_2[i] = cuda_integral_2(cuda_U,h)
        end
    end
    return(
        x,
        t,
        capture_times_flag ? U_set : Array(cuda_U),
        filtration_flag ? (cumsum(I1_dissipated), cumsum(I2_dissipated)) : (nothing, nothing),
        pulse_maximum_flag ? pulse_maximum : nothing,
        integrals_flag ? (I_1, I_2) : (nothing, nothing),
    )
end
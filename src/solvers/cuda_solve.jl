function prepare_filtration_args(filtration_args, N_t)
    @unpack filtration_flag, filtration_tau, filtration_factor, filtration_end_t, filtration_l_nominal = filtration_args
    return(
        filtration_flag ? filtration_tau : nothing,
        filtration_flag ? filtration_factor : nothing,
        filtration_flag ? filtration_end_t : nothing,
        filtration_flag ? [filtration_tau] : nothing,
        filtration_flag ? filtration_l_nominal : nothing,
        filtration_flag ? zeros(N_t) : nothing,
        filtration_flag ? zeros(N_t) : nothing,
    )    
end
function prepare_capture_times_args(capture_times, tspan, N)
    (any(capture_times.>tspan[2]) || any(capture_times.<tspan[1])) && @error "capture_times contains an element out of modeling time"
    capture_times_flag=~isempty(capture_times)
    return(
        capture_times_flag,
        capture_times_flag ? [Vector{ComplexF64}(undef, N) for _ in 1:length(capture_times)] : nothing,
        capture_times_flag ? [1] : nothing,
        capture_times_flag ? [capture_times..., Inf] : nothing,
        capture_times_flag ? [capture_times[1]] : nothing,
    )   
end
function cuda_cycle_iteration!(
    i::Int,
    h::Real,
    tau::Real,
    cuda_M::CuArray{ComplexF64, 2, CUDA.Mem.DeviceBuffer},
    cuda_U::CuArray{ComplexF64, 1, CUDA.Mem.DeviceBuffer},
    ε_2::Real,
    ε_3::Real,
    # Return solution at times given
    capture_times_flag::Bool,      
    U_set::Union{Vector{Vector{ComplexF64}}, Nothing},
    capture_times_index_slider::Union{Vector{Int}, Nothing},
    capture_times::Union{Vector{<:Real}, Nothing},
    t_capture_entry::Union{Vector{<:Real}, Nothing},
    # Return integrals
    integrals_flag::Bool,          
    I_1::Union{Vector{<:Real}, Nothing},
    I_2::Union{Vector{<:Real}, Nothing},
    # Return pulse maximum
    pulse_maximum_flag::Bool,      
    pulse_maximum::Union{Vector{<:Real}, Nothing},
    filtration_flag::Union{Bool, Nothing},
    filtration_tau::Union{Real, Nothing},
    filtration_factor::Union{Real, Function, Nothing},
    filtration_end_t::Union{Real, Nothing},
    filtration_t_slider::Union{Vector{<:Real}, Nothing},
    filtration_l_nominal::Union{Real, Nothing},
    I1_dissipated::Union{Vector{<:Real}, Nothing},
    I2_dissipated::Union{Vector{<:Real}, Nothing},
    # Live plot
    live_plot_solution_flag::Bool,
    live_plot_tau::Real,
    live_plot_t_slider::Union{Vector{<:Real}, Nothing},
    live_t_observable::Union{Nothing, Observable{Float64}},
    observable_U::Union{Nothing, Observable{T}},
    io,
) where T <: Vector{<:Real}
    t_current=(i-1)*tau
    if live_plot_solution_flag && (t_current ≥ live_plot_t_slider[1])
        observable_U[] = abs.(Array(cuda_U))
        live_t_observable[] = t_current
        recordframe!(io)
        live_plot_t_slider[1] += live_plot_tau
    end
    if capture_times_flag && t_capture_entry[1]≠Inf
        if t_current ≥ t_capture_entry[1]
            U_set[capture_times_index_slider[1]] = Array(cuda_U) #Не копировать!
            capture_times_index_slider[1]+=1
            t_capture_entry[1]=capture_times[capture_times_index_slider[1]]
        end
    end
    if pulse_maximum_flag
        pulse_maximum[i]=cuda_maximum(cuda_U)
    end
    if filtration_flag && (t_current <= filtration_end_t)
        if t_current ≥ filtration_t_slider[1]
            filtration_t_slider[1]+=filtration_tau
            cuda_U, (power_I1, power_I2) = cuda_filtration(
                cuda_U,
                h,
                isa(filtration_factor, Function) ? filtration_factor(t_current) : filtration_factor,
                filtration_l_nominal,
            )
            I1_dissipated[i] = power_I1
            I2_dissipated[i] = power_I2
        end
    end
    if integrals_flag
        I_1[i] = cuda_integral_1(cuda_U,h)
        I_2[i] = cuda_integral_2(cuda_U,h)
    end

    #cuda_U[1:end] = cuda_matrix_vector_multiplication(cuda_M,cuda_calculate_V(cuda_U, tau, ε_2, ε_3))
    cuda_U[1:end] = cuda_calculate_V_and_multiplicate(cuda_M, cuda_U, tau, ε_2, ε_3)
end
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
    filtration_args::NamedTuple = (
        filtration_flag = false,
        filtration_tau = 10.0,
        filtration_factor = 1.0,
        filtration_end_t = tspan[2],
        filtration_l_nominal = 100.0,
    ),
    # pulse_maximum calculations
    pulse_maximum_flag = false,
    # record integrals
    integrals_flag = false,
    # times of interest
    capture_times = [],
    # live plot
    live_plot_solution_flag = false,
    live_plot_tau = 1.0,
    # no progress bar
    omit_progress_flag = false,
)
    @unpack filtration_flag, filtration_tau, filtration_factor, filtration_end_t, filtration_l_nominal = filtration_args
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
        @warn "Initial pulse out of x interval."
    
    ε_tail = 1e-5
    (abs(U[1]) < ε_tail) & (abs(U[end]) < ε_tail)  ||
    @warn "Pulse tail modulus exceeds $ε_tail with the value of $(max(abs(U[1]), abs(U[end]))). Consider larger x-interval."

    cuda_U = CuArray(U)
    cuda_M = CuArray(M)

    filtration_tau, filtration_factor, filtration_end_t, filtration_t_slider,
        filtration_l_nominal, I1_dissipated, I2_dissipated = prepare_filtration_args(filtration_args, N_t)

    capture_times_flag, U_set, capture_times_index_slider, capture_times,
        t_capture_entry = prepare_capture_times_args(capture_times, tspan, N)

    pulse_maximum =  pulse_maximum_flag ? zeros(N_t) : nothing

    I_1 = I_2 = integrals_flag ? zeros(N_t) : nothing

    if live_plot_solution_flag
        observable_U = Observable(abs.(U))
        observable_x = Observable(x)
        live_plot_t_slider = [0.0]
        live_t_observable=Observable(live_plot_t_slider[1])
        live_plot = Figure()
        live_axis = Axis(live_plot[1, 1], limits = (minimum(x), maximum(x), 0.0, 1.0))
        lines!(live_axis, observable_x, observable_U)
        #annotate!([(100.0,0.5, Plots.text("t=$(live_t_observable)"))]) TODO
        display(live_plot)
    end

    iteration_args = (
        h,
        tau,
        cuda_M,
        cuda_U,
        ε_2,
        ε_3,
        capture_times_flag, U_set, capture_times_index_slider, capture_times, t_capture_entry,
        integrals_flag, I_1, I_2,
        pulse_maximum_flag, pulse_maximum,
        filtration_flag, filtration_tau, filtration_factor, filtration_end_t, filtration_t_slider,
            filtration_l_nominal, I1_dissipated, I2_dissipated,
        live_plot_solution_flag,
        live_plot_tau,
        live_plot_solution_flag ? live_plot_t_slider : [0.0],
        live_plot_solution_flag ? live_t_observable : nothing,
        live_plot_solution_flag ? observable_U : nothing,
    )

    if live_plot_solution_flag 
        GLMakie.record(live_plot, "output.mp4", framerate=30) do io
            if omit_progress_flag
                for i in 1:N_t
                    cuda_cycle_iteration!(i, iteration_args..., io)
                end
            else
                @showprogress for i in 1:N_t
                    cuda_cycle_iteration!(i, iteration_args..., io)
                end
            end
        end
    else
        io=nothing
        if omit_progress_flag
            for i in 1:N_t
                cuda_cycle_iteration!(i, iteration_args..., io)
            end
        else
            @showprogress for i in 1:N_t
                cuda_cycle_iteration!(i, iteration_args..., io)
            end
        end
    end

    return(
        x,
        t,
        capture_times_flag ? U_set : Array(cuda_U),
        filtration_flag ? (cumsum(I1_dissipated), cumsum(I2_dissipated)) : (zeros(N_t), zeros(N_t)),
        pulse_maximum_flag ? pulse_maximum : nothing,
        integrals_flag ? (I_1, I_2) : (nothing, nothing),
    )
end
# include("run/run_4_fifth_nonlinear_influence.jl")
using Revise, LaTeXStrings
using GeneralizedSchrodingerResearch.BaseUtilities: smooth_vector, relative_error_to_amplitude
using GeneralizedSchrodingerResearch.Solvers: solve, cuda_solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution, construct_approximate_NSE_3_5_7_solution
using Plots, Statistics

k = 0.15
ω = 0.4

ε_2 = 0.5
ε_3 = 0.0

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0

h=0.2
tau=h^2
tspan = (0.0,50000.0)
xspan = (-80.0, 80.0)
capture_times = collect(0.0:1.0:tspan[2])

initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)

filtration_factor, filtration_end_t = 1.02, Inf#tspan[2]
filtration_factor_function = t -> t > filtration_end_t ? 1.0 : ((1.0 - filtration_factor) / filtration_end_t) * t + filtration_factor

x, t, U_set, (I1_dissipated, I2_dissipated), pulse_maximum, (I1, I2) = cuda_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function;
    method = "fourier",
    ε_2 = ε_2,
    ε_3 = ε_3,
    filtration_flag = true,
    filtration_time = 0.8,
    filtration_factor = filtration_factor_function,
    filtration_end_t = filtration_end_t,
    l_nominal = 60.0,
    pulse_maximum_flag = true,
    integrals_flag = true,
    capture_times = capture_times,
    live_plot_solution_flag = false,
    live_plot_tau = 1.0,
); nothing

# Постпроцессинг
begin

    # Вычисление ошибок
    begin
        true_tolerance = Vector{Float64}([])
        possible_NSE_3_5_solution = []
        for i in eachindex(U_set)
            push!(
                possible_NSE_3_5_solution,
                construct_approximate_NSE_3_5_solution(
                    x,
                    h,
                    U_set[i],
                    ε_2,
                    xspan[2]-xspan[1];
                ),
            )
            push!(true_tolerance, relative_error_to_amplitude(possible_NSE_3_5_solution[i],U_set[i]))
        end
        println("В конечный момент ошибка между решениями составляет $(round(true_tolerance[end],digits=5)) %")
    end

    # Построение графика ошибки
    begin
        N_smooth=200
        N_smooth>length(true_tolerance) ? (@warn "Сглаживание true_tolerance задано некорректно и применяться не будет") : 
        (
            n_smooth_from=1; # cld(length(true_tolerance),2);
            true_tolerance[n_smooth_from:end]=smooth_vector(true_tolerance[n_smooth_from:end],N_smooth)
        )
        plot_0 = plot(
            capture_times,
            true_tolerance,
            label="относительная ошибка",
            legendfontsize=11,
            line=(:path,:solid,:black,1),
            dpi=600,
            xlabel="t",
            ylabel=L"\Delta_{\%}^{n}",
            legend=:topleft,
            xlabelfontsize=14, 
            ylabelfontsize=16,
        )
        plot!(
            twinx(),
            capture_times,
            (filtration_factor_function.(capture_times)),
            fillrange=0,
            alpha=0.1,
            fillcolor=:red,
            label=L"k_{f}\left(t\right)",
            lw=1,
            color=:red,
            right_axis=true,
            legend=:topright,
            ylims=(1, 1.03),
            ylabel=L"k_{f}\left(t\right)",
            ylabelfontsize=16,
            legendfontsize=11,
        )
        plot!(
            [t[1], t[end]],
            [0.0, 0.0];
            label="",
            line=(:path,:dash,:black,1), 
        )
        if filtration_end_t < tspan[2]
            plot!(
                [filtration_end_t],
                [0.0];
                label="",
                marker=(:xcross,4,:red,:red)
            )
            annotate!([(filtration_end_t, 1.5, text("отключение\nфильтрации",9))])
        end
        xaxis = Plots.get_axis(Plots.get_subplot(plot_0,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_0,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        plot!(plot_0, tickfontsize=10, yguidefontrotation=0.0)
        savefig(plot_0, "run/plots/run_4/true_tolerance.png")
    end

    # Построение профилей аналитического, численного и начального решений в интересующие моменты времени
    for i in [10,length(U_set)] #eachindex(U_set)
        n_reduce=Int(round(0.00*length(x)))
        n_circshift=500
        plot_1 = plot(
            x[n_reduce+1:end-n_reduce],
            circshift(abs.(U_set[i])[n_reduce+1:end-n_reduce],n_circshift);
            label="численное решение",
            line=(:path,:solid,:blue,2),
            dpi=600,
        )
        plot!(
            x[n_reduce+1:end-n_reduce],
            circshift(
                100 .* (abs.(U_set[i])[n_reduce+1:end-n_reduce] .- abs.(possible_NSE_3_5_solution[i])[n_reduce+1:end-n_reduce]),
                n_circshift
            );
            label="абсолютная ошибка\nмежду решениями * C",
            line=(:path,:solid,:gray,1),
            dpi=600,
        )
        plot!(
            x[n_reduce+1:end-n_reduce],
            abs.(initial_function.(x))[n_reduce+1:end-n_reduce];
            label="начальный импульс",
            line=(:path,:dash,:red,2),
            dpi=600,
            #marker=(:circle,4,:black,:black)
        )
        plot!(
            x[n_reduce+1:end-n_reduce],
            circshift(abs.(possible_NSE_3_5_solution[i])[n_reduce+1:end-n_reduce],n_circshift);
            label="соответствующее\nаналитическое решение",
            line=(:path,:dashdot,:green,3),
            dpi=600,
        )
        xaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        plot!(legend=:topright, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
        xlabel!("x")
        ylabel!("|U|")
        savefig(plot_1, "run/plots/run_4/profiles_T=$(capture_times[i]).png")
    end

    # Построение первого интеграла
    begin
        err_1=(maximum(I1+I1_dissipated) - minimum(I1+I1_dissipated)) / mean(I1+I1_dissipated) * 100
        println("Рассеяно ", (I1[1] - I1[end])/I1[1] * 100, "% от I_1")
        plot_2 = plot(
            t,
            I1;
            dpi=800,
            tickfontsize=10,
            legendfontsize=12,
            guidefontsize=14,
            line=(:path,:solid,:black,2), 
            yformatter = y -> string(round(y,digits=13)),
            label=L"I_{1}",
            legend=:right,
        )
        annotate!(
            [(
                (tspan[2]-tspan[1])/2 + tspan[1],
                (maximum(I1) - minimum(I1))*0.9+minimum(I1),
                text("Рассеяно $(
                    round(
                        (I1[1] - I1[end])/I1[1] * 100,digits=5)
                ) % от " * L"I_{1}",9)
            )]
        )
        xlabel!("t")
        ylabel!("I₁, δ = $(round(err_1,digits=12)) %")
        xaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        if I1_dissipated != nothing
            plot!(
                t,
                I1+I1_dissipated;
                tickfontsize=10,
                legendfontsize=12,
                guidefontsize=14,
                label=L"I_{1} + I_{1f}",
                line=(:path,:solid,:green,2), 
            )
        end
        savefig(plot_2, "run/plots/run_4/I1_influenced.png")
    end

    # Построение второго интеграла
    begin
        err_2=(maximum(I2+I2_dissipated) - minimum(I2+I2_dissipated)) / mean(I2+I2_dissipated) * 100
        println("Рассеяно ", (I2[1] - I2[end])/I2[1] * 100, "% от I_2")
        plot_3 = plot(
            t,
            smooth_vector(I2,100); 
            dpi=800,
            tickfontsize=10,
            legendfontsize=12,
            guidefontsize=14,
            line=(:path,:solid,:black,2), 
            yformatter = y -> string(round(y,digits=5)),
            label=L"I_{2}",
            legend=:right,
        )
        annotate!(
            [(
                (tspan[2]-tspan[1])/2 + tspan[1],
                (maximum(I2) - minimum(I2))*0.9+minimum(I2),
                text("Рассеяно $(round((I2[1] - I2[end])/I2[1] * 100,digits=5)) % от "*L"I_{2}",9)
            )]
        )
        xlabel!("t")
        ylabel!("I₂, δ = $(round(err_2,digits=6)) %")
        xaxis = Plots.get_axis(Plots.get_subplot(plot_3,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_3,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        if I2_dissipated != nothing
            plot!(
                t,
                smooth_vector(I2+I2_dissipated,100);
                tickfontsize=10,
                legendfontsize=12,
                guidefontsize=14,
                label=L"I_{2} + I_{2f}",
                line=(:path,:solid,:green,2), 
            )
        end
        savefig(plot_3, "run/plots/run_4/I2_influenced.png")
    end

    # Построение поведения пика импульса
    if !isnothing(pulse_maximum)
        plot_4 = plot(
            t,
            smooth_vector(pulse_maximum,600); 
            line=(:path,:solid,:black,1), 
            legend = false,
            dpi=800,
            tickfontsize=10,
            guidefontsize=14,
        )
        xaxis = Plots.get_axis(Plots.get_subplot(plot_4,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_4,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        xlabel!("t")
        ylabel!(L"\max|U^{m}|")
        savefig(plot_4, "run/plots/run_4/pike_behaviour.png")

        time_gap=100.0
        t1_index=Int(round(time_gap/tau))
        t2_index=Int(round((tspan[2]-time_gap)/tau))
        plot_5 = plot(
            t[1:t1_index],
            pulse_maximum[1:t1_index]; 
            line=(:path,:solid,:black,1), 
            legend = false,
            dpi=800,
            tickfontsize=10,
            guidefontsize=14,
            margin=2Plots.mm,
        )
        xaxis = Plots.get_axis(Plots.get_subplot(plot_5,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_5,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        xlabel!("t")
        ylabel!(L"\max|U^{m}|")
        savefig(plot_5, "run/plots/run_4/pike_behaviour_initial.png")

        N_smooth=200
        plot_6 = plot(
            t[t2_index:(end-N_smooth)],
            smooth_vector(pulse_maximum[t2_index:end],N_smooth)[1:(end-N_smooth)]; 
            line=(:path,:solid,:black,1), 
            legend = false,
            dpi=800,
            tickfontsize=10,
            guidefontsize=14,
            margin=2Plots.mm,
        )
        xaxis = Plots.get_axis(Plots.get_subplot(plot_6,1),:x)
        yaxis = Plots.get_axis(Plots.get_subplot(plot_6,1),:y)
        xaxis[:gridalpha] = 0.4
        yaxis[:gridalpha] = 0.4
        xlabel!("t")
        ylabel!(L"\max|U^{m}|")
        savefig(plot_6, "run/plots/run_4/pike_behaviour_final.png")
    end

    # Сохранение результатов
    #= begin
        U_set_string=join(U_set, ", ")
        x_string=join(x, ", ")
        t_string = join(t, ", ")
        I1_dissipated_string = join(I1_dissipated, ", ")
        I2_dissipated_string = join(I2_dissipated, ", ")
        pulse_maximum_string = join(pulse_maximum, ", ")
        I1_string = join(I1, ", ")
        I2_string = join(I2, ", ")
        if !isdir("run/results")
            mkdir("run/results")
        end
        open("run/results/x.txt", "w") do file
            write(file, x_string)
        end
        open("run/results/U_set.txt", "w") do file
            write(file, U_set_string)
        end
        open("run/results/t.txt", "w") do file
            write(file, t_string)
        end
        open("run/results/I1_dissipated.txt", "w") do file
            write(file, I1_dissipated_string)
        end
        open("run/results/I2_dissipated.txt", "w") do file
            write(file, I2_dissipated_string)
        end
        open("run/results/pulse_maximum.txt", "w") do file
            write(file, pulse_maximum_string)
        end
        open("run/results/I1.txt", "w") do file
            write(file, I1_string)
        end
        open("run/results/I2.txt", "w") do file
            write(file, I2_string)
        end
    end =#

end
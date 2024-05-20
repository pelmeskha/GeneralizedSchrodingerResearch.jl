# include("run/run_4_nonlinear_influence.jl")
using Revise, LaTeXStrings
using GeneralizedSchrodingerResearch.BaseUtilities: smooth_vector
using GeneralizedSchrodingerResearch.Solvers: solve, cuda_solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution, construct_approximate_NSE_3_5_7_solution
using Plots, Statistics

k = 0.15
ω = 0.4

ε_2 = 0.65
ε_3 = 0.08

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0

h=0.25
tau=h^2
tspan = (0.0,8000.0)
xspan = (-80.0, 80.0)
capture_times = [10.0, tspan[2]]

initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)

filtration_factor, filtration_end_t = 1.03, 8000.0
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
    live_plot_solution_flag = true,
    live_plot_tau = 1.0,
)

naive_tolerance = (pulse_maximum.-maximum(abs.(U_set[end]))) ./ maximum(abs.(U_set[end])) .* 100
possible_NSE_3_5_solution = []
for i in eachindex(U_set)
    push!(
        possible_NSE_3_5_solution,
        construct_approximate_NSE_3_5_solution(
            x,
            capture_times[i],
            U_set[i],
            k,
            ε_2,
            theta_0,
            z_0,
            xspan[2]-xspan[1],
        ),
    )
end

possible_NSE_3_5_7_solution = []
for i in eachindex(U_set)
    push!(
        possible_NSE_3_5_7_solution,
        construct_approximate_NSE_3_5_7_solution(
            x,
            capture_times[i],
            U_set[i],
            k,
            ε_2,
            theta_0,
            z_0,
            ξ₀,
            xspan[2]-xspan[1],
        ),
    )
end


for i in eachindex(U_set)
    reduce=Int(round(0.00*length(x)))+1
    plot_1 = plot(
        x[reduce:end-reduce],
        abs.(U_set[i])[reduce:end-reduce];
        label="численное решение",
        line=(:path,:solid,:blue,2),
        dpi=600,
    )
    plot!(
        x[reduce:end-reduce],
        abs.(initial_function.(x))[reduce:end-reduce];
        label="начальный импульс",
        line=(:path,:dash,:red,2),
        dpi=600,
        #marker=(:circle,4,:black,:black)
    )
    #= plot!(
        x[reduce:end-reduce],
        abs.(possible_NSE_3_5_solution[i])[reduce:end-reduce];
        label="соответствующее\nаналитическое решение",
        lw=3,
        ls=:dashdot,
        dpi=600,
    ) =#
    plot!(
        x[reduce:end-reduce],
        abs.(possible_NSE_3_5_7_solution[i])[reduce:end-reduce];
        label="соответствующее\nаналитическое решение",
        lw=3,
        ls=:dashdot,
        dpi=600,
    )
    xaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:x)
    yaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:y)
    xaxis[:gridalpha] = 0.4
    yaxis[:gridalpha] = 0.4
    plot!(legend=:topright, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
    xlabel!("x")
    ylabel!("|U|")
    savefig(plot_1, "run/plots/profiles_T=$(capture_times[i]).png")
end

err_1=(maximum(I1+I1_dissipated) - minimum(I1+I1_dissipated)) / mean(I1+I1_dissipated) * 100
println("Рассеяно ", ((I1+I1_dissipated)[end] - I1[end])/(I1+I1_dissipated)[end] * 100, "% от I_1")
plot_2 = plot(
    t,
    I1;
    dpi=800,
    tickfontsize=10,
    legendfontsize=10,
    guidefontsize=14,
    line=(:path,:solid,:black,1), 
    yformatter = y -> string(round(y,digits=13)),
    label=L"I_{1}",
    legend=:right,
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
        legendfontsize=10,
        guidefontsize=14,
        label=L"I_{1} + I_{1f}",
        line=(:path,:solid,:green,2), 
    )
end
savefig(plot_2, "run/plots/I1_influenced.png")

err_2=(maximum(I2+I2_dissipated) - minimum(I2+I2_dissipated)) / mean(I2+I2_dissipated) * 100
println("Рассеяно ", ((I2+I2_dissipated)[end] - I2[end])/(I2+I2_dissipated)[end] * 100, "% от I_2")
plot_3 = plot(
    t,
    I2; 
    dpi=800,
    tickfontsize=10,
    legendfontsize=10,
    guidefontsize=14,
    line=(:path,:solid,:black,1), 
    yformatter = y -> string(round(y,digits=5)),
    label=L"I_{2}",
    legend=:right,
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
        I2+I2_dissipated;
        tickfontsize=10,
        legendfontsize=10,
        guidefontsize=14,
        label=L"I_{2} + I_{2f}",
        line=(:path,:solid,:green,2), 
    )
end
savefig(plot_3, "run/plots/I2_influenced.png")

if !isnothing(pulse_maximum)
    plot_4 = plot(
        t,
        smooth_vector(naive_tolerance,200); 
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
    ylabel!(L"\Delta_{\%}^{n}")
    savefig(plot_4, "run/plots/pike_tolerance.png")
end
possible=join(abs.(possible_NSE_3_5_7_solution[end]), " ")
exact=join(abs.(U_set[end]), " ")
x_str=join(x, " ")
if !isdir("run/results")
    mkdir("run/results")
end
open("run/results/x.txt", "w") do file
    write(file, x_str)
end
open("run/results/possible.txt", "w") do file
    write(file, possible)
end
open("run/results/exact.txt", "w") do file
    write(file, exact)
end
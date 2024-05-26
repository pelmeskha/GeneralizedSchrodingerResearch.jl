#include("run/run_6_dissipation_of_epsilon.jl")
using Revise, LaTeXStrings, ProgressMeter
using GeneralizedSchrodingerResearch.BaseUtilities: smooth_vector, relative_error_to_amplitude
using GeneralizedSchrodingerResearch.Solvers: solve, cuda_solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution, construct_approximate_NSE_3_5_7_solution
using Plots, Statistics

k = 0.15
ω = 0.4

#ε_2_vector = collect(-1.0:0.05:0.0)
ε_2_vector = collect(0.0:0.05:1.0)
ε_3 = 0.0

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0
h=0.2
tau=h^2
tspan = (0.0,1000.0)
xspan = (-80.0, 80.0)

initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)

filtration_factor, filtration_end_t = 1.02, tspan[2]
filtration_factor_function = t -> t > filtration_end_t ? 1.0 : ((1.0 - filtration_factor) / filtration_end_t) * t + filtration_factor

I1_dissipated_percents=Vector{Float64}([])
I2_dissipated_percents=Vector{Float64}([])
final_pulse_amplitude=Vector{Float64}([])
@showprogress for ε_2 in ε_2_vector
    _, t, _, (I1_dissipated, I2_dissipated), pulse_maximum, (I1, I2) = cuda_solve(
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
        live_plot_solution_flag = false,
        live_plot_tau = 5.0,
    )
    push!(I1_dissipated_percents, ((I1+I1_dissipated)[end] - I1[end])/(I1+I1_dissipated)[end] * 100)
    push!(I2_dissipated_percents, ((I2+I2_dissipated)[end] - I2[end])/(I2+I2_dissipated)[end] * 100)
    push!(final_pulse_amplitude, pulse_maximum[end])
end

plot_1 = plot(
    ε_2_vector,
    I1_dissipated_percents;
    dpi=800,
    tickfontsize=10,
    legendfontsize=12,
    guidefontsize=14,
    line=(:path,:solid,:red,1), 
    yformatter = y -> string(round(y,digits=13)),
    label=L"\Delta I_{1} \%",
    marker=(:utriangle,3,:red,:red),
    legend=:topright,
    yticks=7,
)
plot!(
    ε_2_vector,
    I2_dissipated_percents;
    line=(:path,:solid,:blue,1), 
    yformatter = y -> string(round(y,digits=13)),
    label=L"\Delta I_{2} \%",
    marker=(:dtriangle,3,:blue,:blue),
)
xlabel!(L"\varepsilon_{2}")
ylabel!("Потери из-за нелинейности, %")
xaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:y)
xaxis[:gridalpha] = 0.4
yaxis[:gridalpha] = 0.4
savefig(plot_1, "run/plots/dissipation_over_epsilon_2.png")

plot_2 = plot(
    ε_2_vector,
    final_pulse_amplitude;
    dpi=800,
    tickfontsize=10,
    legendfontsize=12,
    guidefontsize=14,
    line=(:path,:solid,:black,1), 
    yformatter = y -> string(round(y,digits=13)),
    label="",
    marker=(:circle,1,:black,:black),
)
xlabel!(L"\varepsilon_{2}")
ylabel!(L"\max|U_{N_{t}}|")
xaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:y)
xaxis[:gridalpha] = 0.4
yaxis[:gridalpha] = 0.4
savefig(plot_2, "run/plots/final_amplitude_over_epsilon_2.png")

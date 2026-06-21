#include("run/run_8_nonlinear_diagram.jl")
using Revise, LaTeXStrings, ProgressMeter, DelimitedFiles, Plots, Statistics, CSV, DataFrames
using GeneralizedSchrodingerResearch.BaseUtilities: smooth_vector, relative_error_to_amplitude
using GeneralizedSchrodingerResearch.Solvers: solve, cuda_solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution, 
    construct_approximate_NSE_3_5_7_solution, construct_best_approximate_NSE_3_5_7_solution

k = 0.15
ω = 0.6

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0

h=0.2
tau=h^2
tspan = (0.0,100.0)
xspan = (-60.0, 60.0)

initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)

filtration_factor, filtration_end_t = 1.00, 0.0
filtration_factor_function = t -> t > filtration_end_t ? 1.0 : ((1.0 - filtration_factor) / filtration_end_t) * t + filtration_factor

filtration_args=(
    filtration_flag=false,
    filtration_tau=0.8,
    filtration_factor=filtration_factor_function,
    filtration_end_t=filtration_end_t,
    filtration_l_nominal=60.0,
)

ε_2_vector = collect(-1.0:0.1:1.0)
ε_3_vector = collect(-0.5:0.05:0.5)
average_tolerance = zeros(length(ε_2_vector), length(ε_3_vector))

p_outer = Progress(length(ε_2_vector), desc="Outer loop progress")
p_inner = Progress(length(ε_3_vector), desc="Inner loop progress")

total_iterations = length(ε_2_vector) * length(ε_3_vector)
p = Progress(total_iterations)

Base.Threads.nthreads() = 4

@showprogress for index_ε_2 in eachindex(ε_2_vector)
    for index_ε_3 in eachindex(ε_3_vector)
        _, _, _, (I1_dissipated, I2_dissipated), _, (I1, I2) = cuda_solve(
            tspan,
            xspan,
            tau,
            h,
            initial_function;
            method = "fourier",
            ε_2 = ε_2_vector[index_ε_2],
            ε_3 = ε_3_vector[index_ε_3],
            filtration_args,
            pulse_maximum_flag = false,
            integrals_flag = true,
        )
        local err_1=abs((maximum(I1+I1_dissipated) - minimum(I1+I1_dissipated)) / mean(I1+I1_dissipated) * 100)
        local err_2=abs((maximum(I2+I2_dissipated) - minimum(I2+I2_dissipated)) / mean(I2+I2_dissipated) * 100)
        average_tolerance[index_ε_2, index_ε_3] = maximum([err_1, err_2])
        next!(p)
    end
end

finish!(p)

average_tolerance_ = map(x -> x > 100.0 ? 100.0 : x, average_tolerance)
heatmap(
    ε_2_vector,
    ε_3_vector,
    average_tolerance_,
    #color=:gray,
    xlabel=L"\varepsilon_{2}",
    ylabel=L"\varepsilon_{3}",
    title="Нарушение законов сохранения, %",
    legend=false,
    tickfontsize=10,
    legendfontsize=12,
    guidefontsize=14,
)
vline!([0], linestyle=:dash, color=:red)
hline!([0], linestyle=:dash, color=:blue)
plot!(;dpi=600)
savefig("run/plots/run_8/diagram(no_filt), T=$(tspan[2]), h=$h, k=$k, ω=$ω.png")

# Сохранение результатов
#= begin
    ε_2_vector_string=join(ε_2_vector, ", ")
    ε_3_vector_string=join(ε_3_vector, ", ")
    average_tolerance_string = join(average_tolerance, ", ")
    if !isdir("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω")
        mkdir("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω")
    end
    open("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω/ε_2_vector.txt", "w") do file
        write(file, ε_2_vector_string)
    end
    open("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω/ε_3_vector.txt", "w") do file
        write(file, ε_3_vector_string)
    end
    open("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω/average_tolerance.txt", "w") do file
        write(file, average_tolerance_string)
    end
    writedlm("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω/average_tolerance_2.txt", average_tolerance)
    average_tolerance_df = DataFrame(average_tolerance, :auto)
    CSV.write("run/results/diagram_T=$(tspan[2]), h=$h, k=$k, ω=$ω/average_tolerance.csv", average_tolerance_df)
end =#
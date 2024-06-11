#include("run/run_6_postprocess_numerical_solution.jl")
using Revise, LaTeXStrings, Plots, Statistics, Interpolations, Optim
using GeneralizedSchrodingerResearch.BaseUtilities: smooth_vector,
    parse_to_vector, set_grid_alpha, shift_pulse_to_center, find_threshold
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution, 
    construct_approximate_NSE_3_5_7_solution

numerical_solution_filename = "run/results/exact.txt"
x_nodes_filename = "run/results/x.txt"
numerical_solution = parse_to_vector(numerical_solution_filename)
x_nodes = parse_to_vector(x_nodes_filename)

abs_U = abs.(numerical_solution)
bottom_gap = find_threshold(abs_U , 12)
x, y = shift_pulse_to_center(x_nodes, abs_U ; y_threshold=maximum(abs_U )-bottom_gap)

h=round(x[2]-x[1], digits=12) # Костыльно
x_range=x[1]:h:x[end]
interpolator = cubic_spline_interpolation((x_range), y)
interpolator_function = xi -> interpolator(xi)
minimization_function = xi -> -interpolator(xi)

x_subgrid = collect(range(x[1], stop=x[end], length=100))
y_subgrid = interpolator_function.(x_subgrid)

result = optimize(minimization_function, minimum(x), maximum(x))
x_max, y_max = result.minimizer, -result.minimum
println("Maximum value: $y_max")
#= plot_1 = plot(
    x,
    y;
    line=(:path,:solid,:red,2),
    dpi=600,
    marker=(:circle,1,:black,:black)
)
scatter!(
    [x_max],
    [y_max];
    marker=(:circle,3,:blue,:blue)
)
plot!(
    x_subgrid,
    y_subgrid;
    line=(:path,:dash,:black,1),
) =#

y_possible = construct_approximate_NSE_3_5_7_solution(
    x_nodes,
    h,
    numerical_solution,
    0.65,
    x_nodes[end]-x_nodes[1]+h;
    use_M_values=true,
    M0=-5.4381083,
)

N_shift=801
idxs=664:936
println("средняя ошибка по импульсу: ", mean((circshift(abs.(y_possible),N_shift)[idxs] - circshift(abs.(numerical_solution),N_shift)[idxs]) ./ 
    circshift(abs.(numerical_solution),N_shift)[idxs] .* 100.0), " %")

N=790
idxs=(1+N):(length(x_nodes)-N)

plot_1 = plot(
    x_nodes[idxs],
    circshift(abs.(y_possible),N_shift)[idxs];
    line=(:path,:solid,:blue,1),
    marker=(:circle,1,:red,:red),
    dpi=600,
)
plot!(
    x_nodes[idxs],
    circshift(abs.(numerical_solution),N_shift)[idxs];
    line=(:path,:solid,:red,1),
    marker=(:circle,1,:blue,:blue),
)
set_grid_alpha(plot_1, 0.4)
display(plot_1)
# include("run/run_0_plot_solitons.jl")
using Revise, LaTeXStrings
using GeneralizedSchrodingerResearch.Solvers: solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution
using GeneralizedSchrodingerResearch.BaseUtilities: M0_M1_to_ε2_ε3, ε2_ε3_to_M0_M1
using Plots, Statistics

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0

k = 0.0
ω = 0.0

L = 160.0
h = 0.1
N = Int(L / h)
N_x = N + 1
j = range(-N / 2, stop = N / 2 , length = N_x)
x = collect(j .* h)[1:end-1]


#= initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)
analytical_solution = (x, t) -> NSE_3_soliton(x, t, k, ω, theta_0, z_0) =#

#initial_function = (x) ->NSE_3_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
#analytical_solution = (x, t) -> NSE_3_5_soliton(x, t, k, ω, -1, theta_0, z_0)

ε_2, ε_3 =  0.65, 0.08 #M0_M1_to_ε2_ε3(-4.0, 20.0)
println(ε_2," ", ε_3)
interpolator = precompile_NSE_3_5_7_soliton(ε_2,ε_3,z_0,ξ₀,L)
soliton_1 = (x,t) -> NSE_3_5_7_soliton(x,t,k,ω,theta_0,z_0,interpolator)

#= ε_2, ε_3 = 0.65, 0.08 #M0_M1_to_ε2_ε3(-26.0, 108.0)
println(ε_2," ", ε_3)
interpolator = precompile_NSE_3_5_7_soliton(ε_2,ε_3,z_0,ξ₀,L)
soliton_2 = (x,t) -> NSE_3_5_7_soliton(x,t,k,ω,theta_0,z_0,interpolator) =#

plot_1 = plot(
    x,
    abs.(soliton_1.(x,0.0));
    label="солитон 1",
    line=(:path,:solid,:red,2),
    dpi=600,
    #marker=(:circle,4,:black,:black)
)
#= plot!(
    x,
    abs.(soliton_2.(x,0.0));
    label="солитон 2",
    line=(:path,:dash,:black,2),
    dpi=600,
) =#
xaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:y)
xaxis[:gridalpha] = 0.4
yaxis[:gridalpha] = 0.4
plot!(legend=:topright, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
xlabel!("x")
ylabel!("|U|")
savefig(plot_1, "run/plots/run_7/solitons.png")

# include("run/run_3_solitary_wave_prop.jl")
using Revise, LaTeXStrings
using GeneralizedSchrodingerResearch.Solvers: solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution
using Plots, Statistics

k = 0.15
ω = 0.45

ε_2 = 1.0
ε_3 = 0.3

theta_0 = 0.0
z_0 = 0.0
ξ₀ = 0.0

h=0.2
tau=h^2
tspan = (0.0,500.0)
xspan = (-60.0, 60.0)

initial_function = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)
analytical_solution = (x, t) -> NSE_3_soliton(x, t, k, ω, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1])

#initial_function = (x) ->NSE_3_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
#analytical_solution = (x, t) -> NSE_3_5_soliton(x, t, k, ω, -1, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1])

#= interpolator = precompile_NSE_3_5_7_soliton(ε_2,ε_3,z_0,ξ₀,xspan[2]-xspan[1])
analytical_solution = (x, t) -> NSE_3_5_7_soliton(x,t,k,ω,theta_0,z_0,interpolator; cycle=true, L=xspan[2]-xspan[1])
initial_function = (x) -> analytical_solution(x,0.0) =#

x, t, U, (I1_dissipated, I2_dissipated), tolerance, (I1, I2) = solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function;
    method = "fourier",
    ε_2 = ε_2,
    ε_3 = ε_3,
    filtration_flag = false,
    filtration_time = 2.0,
    filtration_factor = 1 + 2e-2,
    l_nominal = 60.0,
    tolerance_flag = true,
    analytical_solution = analytical_solution,
    integrals_flag = true,
)

#= =possible_NSE_5_solution = construct_approximate_NSE_3_5_solution(
    x,
    U,
    ε_2,
    theta_0,
    z_0,
) =#

reduce=Int(round(0.00*length(x)))+1
plot_1 = plot(
    x[reduce:end-reduce],
    abs.(U)[reduce:end-reduce];
    label="численное решение",
    line=(:path,:solid,:black,2),
    dpi=800,
)
plot!(
    x[reduce:end-reduce],
    abs.(analytical_solution.(x,0))[reduce:end-reduce];
    label="начальный импульс",
    line=(:path,:solid,:red,2),
    #marker=(:circle,4,:black,:black)
)
plot!(
    x[reduce:end-reduce],
    abs.(analytical_solution.(x,tspan[2]))[reduce:end-reduce];
    label="аналитическое решение",
    line=(:dash,:dash,:green,4),
    #= markershape = :d,
    markersize = 1,
    markercolor = :green,
    markerstrokewidth = 3,
    markerstrokecolor = :green, =#
)
#plot!(x[reduce:end-reduce],abs.(possible_NSE_5_solution)[reduce:end-reduce]; label="соответствующее\nаналитическое решение", lw=3, ls=:dashdot)
plot!(legend=:topright, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
xlabel!("x")
ylabel!("|U|")
savefig(plot_1, "run/plots/run_3/solution_profiles.png")

err_1=(maximum(I1) - minimum(I1)) / mean(I1) * 100
y_ticks= [4.8-2e-12, 4.8-1e-12, 4.8, 4.8+1e-12]
ytick_labels = ["4.8 - 2e-12", "4.8 - 1e-12", "4.8", "4.8 + 1e-12"]
plot_2 = plot(
    t,
    I1;
    legend = false,
    dpi=800,
    tickfontsize=10,
    guidefontsize=14,
    line=(:path,:solid,:black,1), 
    #yformatter = y -> string(round(y,digits=13)),
    #yticks=(y_ticks, ytick_labels)
)
xlabel!("t")
ylabel!("I₁, δ = $(round(err_1,digits=6)) %")
#ylabel!("I₁, δ = $(round(err_1,digits=6)) %")
xaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_2,1),:y)
xaxis[:gridalpha] = 0.4
yaxis[:gridalpha] = 0.4
if I1_dissipated != nothing
    plot!(t,I1+I1_dissipated; label="I1 + I1_dissipated")
end
#plot!(legend=:outerbottom)
savefig(plot_2, "run/plots/run_3/I1.png")

err_2=(maximum(I2) - minimum(I2)) / mean(I2) * 100
plot_3 = plot(
    t,
    I2; 
    legend = false,
    dpi=800,
    tickfontsize=10,
    guidefontsize=14,
    line=(:path,:solid,:black,1), 
    yformatter = y -> string(round(y,digits=5)),
)
xlabel!("t")
ylabel!("I₂, δ = $(round(err_2,digits=6)) %")
xaxis = Plots.get_axis(Plots.get_subplot(plot_3,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_3,1),:y)
xaxis[:gridalpha] = 0.4
yaxis[:gridalpha] = 0.4
if I2_dissipated != nothing
    plot!(t,I2+I2_dissipated; label="I2 + I2_dissipated")
end
savefig(plot_3, "run/plots/run_3/I2.png")

if !isnothing(tolerance)
    plot_4 = plot(
        t,
        abs.(tolerance); 
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
    savefig(plot_4, "run/plots/run_3/tolerance.png")
end
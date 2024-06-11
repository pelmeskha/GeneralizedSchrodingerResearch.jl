# include("run/run_2_step_error.jl")
using Revise, LaTeXStrings
using GeneralizedSchrodingerResearch.Solvers: solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_3_5_soliton,                  
    NSE_3_soliton, precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_3_5_solution
using Plots, Statistics

k = 1.6
ω = 0.4

ε_2 = 2.16
ε_3 = 0.99

theta_0 = 0.0
z_0 = -20.0
ξ₀ = 0.0

tspan = (0.0,16.0)
xspan = (-80.0, 80.0)

h_list=[1.0, 0.8, 0.5, 0.4, 0.25]
tau_list=h_list.^2

Fourier_error=zeros(length(h_list))
Fourier_time=zeros(length(h_list))

interpolator = precompile_NSE_3_5_7_soliton(ε_2,ε_3,z_0,ξ₀,xspan[2]-xspan[1])
analytical_solution = (x, t) -> NSE_3_5_7_soliton(x,t,k,ω,theta_0,z_0,interpolator; cycle=true, L=xspan[2]-xspan[1])
initial_function = (x) -> analytical_solution(x,0.0)

for i in eachindex(h_list)
    local h=h_list[i]
    local tau=tau_list[i]
    _, _, _, (_, _), tolerance, (_, _) = solve(
        tspan,
        xspan,
        tau,
        h,
        initial_function;
        method = "fourier",
        ε_2 = ε_2,
        ε_3 = ε_3,
        tolerance_flag = true,
        analytical_solution = analytical_solution,
        integrals_flag = true,
    )
    Fourier_error[i]=maximum(tolerance)
end

plot_1 = plot(
    (h_list),
    (Fourier_error);
    line=(:path,:solid,:black,2),
    marker=(:utriangle,5,:black,:gray),
    dpi=800,
    tickfontsize=10,
    guidefontsize=14,
    legend = false,
)
xlabel!("h")
ylabel!(L"\Delta_{\%}")
plot!(yguidefontrotation=0.0)
xaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:x)
yaxis = Plots.get_axis(Plots.get_subplot(plot_1,1),:y)
xaxis[:gridalpha] = 0.4
xaxis[:flip] = true
yaxis[:gridalpha] = 0.4
savefig(plot_1, "run/plots/run_2/step_error.png")

# include("run/solitary_wave_prop.jl")
using Revise
using GeneralizedSchrodingerResearch.Solvers: solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_5_soliton, NSE_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_5_solution
using Plots, Statistics

k = 0.15
ω = 0.4

theta_0 = 0.0
z_0 = 0.0

h=0.25
tau=h^2 /10
tspan = (0.0,25.0)
xspan = (-80.0, 80.0)

initial_function_3 = (x) -> NSE_soliton(x, 0.0, k, ω, theta_0, z_0)
analytical_solution_3 = (x, t) -> NSE_soliton(x, t, k, ω, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
#initial_function_5 = (x) ->NSE_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
#analytical_solution_5 = (x, t) -> NSE_5_soliton(x, t, k, ω, -1, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
ε_2 = -0.0
x, t, U, (I1_dissipated, I2_dissipated), tolerance, (I1, I2) = solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function_3;
    method = "fourier",
    ε_2 = ε_2,
    ε_3 = 0.0,
    filtration_flag = false,
    filtration_time = 2.0,
    filtration_factor = 1 + 2e-2,
    l_nominal = 60.0,
    tolerance_flag = true,
    analytical_solution = analytical_solution_3,
    integrals_flag = true,
)

#= possible_NSE_5_solution = construct_approximate_NSE_5_solution(
    x,
    U,
    ε_2,
    theta_0,
    z_0,
) =#

reduce=Int(round(0.00*length(x)))+1
plot_1 = plot(x[reduce:end-reduce],abs.(U)[reduce:end-reduce]; label="численное решение", line=(:path,:solid,:blue,2))
plot!(x[reduce:end-reduce],abs.(analytical_solution_3.(x,0))[reduce:end-reduce]; label="начальный импульс", lw=2, ls=:dot)
#plot!(x[reduce:end-reduce],abs.(possible_NSE_5_solution)[reduce:end-reduce]; label="соответствующее\nаналитическое решение", lw=3, ls=:dashdot)
plot!(legend=:topright, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
xlabel!("x")
ylabel!("|U|")
savefig(plot_1, "run/solution_profiles.png")

err_1=(maximum(I1) - minimum(I1)) / mean(I1) * 100
plot_2 = plot(t,I1; label="First integral, rel_err = $err_1")
if I1_dissipated != nothing
    plot!(t,I1+I1_dissipated; label="I1 + I1_dissipated")
end
plot!(legend=:outerbottom)
savefig(plot_2, "run/I1.png")

err_2=(maximum(I2) - minimum(I2)) / mean(I2) * 100
plot_3 = plot(t,I2; label="Second integral, rel_err = $err_2")
if I2_dissipated != nothing
    plot!(t,I2+I2_dissipated; label="I2 + I2_dissipated")
end
plot!(legend=:outerbottom)
savefig(plot_3, "run/I2.png")

if !isnothing(tolerance)
    plot_4 = plot(t,abs.(tolerance); line=(:path,:solid,:black,1), legend = false)
    xlabel!("t")
    ylabel!("относительная ошибка, %")
    savefig(plot_4, "run/tolerance.png")
end
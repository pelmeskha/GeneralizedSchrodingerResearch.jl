# include("run/solitary_wave_prop.jl")
using Revise
using GeneralizedSchrodingerResearch.Solvers: fourier_solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_5_soliton, NSE_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_5_solution
using Plots, Statistics

k = 0.1
ω = 0.2

theta_0 = 0.0
z_0 = 0.0

tspan = (0.0,1870.0)
xspan = (-80.0, 80.0)
h=0.2
tau=h^2

initial_function_3 = (x) -> NSE_soliton(x, 0.0, k, ω, theta_0, z_0)
analytical_solution_3 = (x, t) -> NSE_soliton(x, t, k, ω, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
#initial_function_5 = (x) ->NSE_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
#analytical_solution_5 = (x, t) -> NSE_5_soliton(x, t, k, ω, -1, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
ε_2 = -0.5
x, t, U, (I1_dissipated, I2_dissipated), tolerance, (I1, I2) = fourier_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function_3;
    ε_2 = ε_2,
    ε_3 = 0.0,
    filtration_flag = true,
    filtration_time = 2.0,
    filtration_factor = 1 + 2e-2,
    l_nominal = 60.0,
    tolerance_flag = false,
    analytical_solution = analytical_solution_3,
    integrals_flag = true,
)

possible_NSE_5_solution = construct_approximate_NSE_5_solution(
    x,
    U,
    ε_2,
    theta_0,
    z_0,
)

plot_1 = plot(x,abs.(U); label="итоговое численное решение")
plot!(x,abs.(analytical_solution_3.(x,0)); label="решение в начальный момент")
plot!(x,abs.(possible_NSE_5_solution); label="подобранное аналитическое решение")
plot!(legend=:outerbottom)
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
# include("run/solitary_wave_prop.jl")
using Revise
using GeneralizedSchrodingerResearch
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_5_soliton, NSE_soliton
using Plots

k = 0.1
ω = 0.2

theta_0 = 0
z_0 = 0

tspan = (0,1870)
xspan = (-80.0, 80.0)
h=0.25
tau=h^2

initial_function_3 = (x) -> NSE_soliton(x, 0, k, ω, theta_0, z_0)
analytical_solution_3 = (x, t) -> NSE_soliton(x, t, k, ω, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
#initial_function_5 = (x) ->NSE_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
#analytical_solution_5 = (x, t) -> NSE_5_soliton(x, t, k, ω, -1, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)
ε_2 = -0.5
x, U, power, tolerance, (I1, I2) = GeneralizedSchrodingerResearch.Solvers.fourier_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function_3;
    ε_2 = ε_2,
    ε_3 = 0,
    filtration_flag = true,
    filtration_step = 2,
    filtration_factor = 1 + 2e-2,
    l_nominal = 60,
    tolerance_flag = false,
    analytical_solution = analytical_solution_3,
    integrals_flag = true,
)

y_max=maximum(abs.(U))
possible_μ = 2/3 * (2 * ε_2 * y_max^4 + 3 * y_max^2)
possible_solution_5 = NSE_5_soliton.(x, 0, 0, possible_μ/4, ε_2, theta_0, z_0)
_z_0 = x[argmax(abs.(possible_solution_5))]-x[argmax(abs.(U))]
possible_solution_5 = NSE_5_soliton.(x, 0, 0, possible_μ/4, ε_2, theta_0, -_z_0)

plot(x,abs.(U); label="итоговое численное решение")
plot!(x,abs.(analytical_solution_3.(x,0)); label="решение в начальный момент")
plot!(x,abs.(possible_solution_5); label="подобранное аналитическое решение")
plot!(legend=:outerbottom)
# include("run/solitary_wave_prop.jl")
using Revise
using GeneralizedSchrodingerResearch
using Plots

k = 0.3
ω = 0.213

theta_0 = 0
z_0 = 0

tspan = (0,1000)
xspan = (-10, 10)
h=0.25
tau=h^2

initial_function = (x) -> GeneralizedSchrodingerResearch.AnalyticalSolutions.NSE_soliton(x, 0, k, ω, theta_0, z_0)
initial_function_5 = (x) -> GeneralizedSchrodingerResearch.AnalyticalSolutions.NSE_5_soliton(x, 0, k, ω, -1, theta_0, z_0)
x, U = GeneralizedSchrodingerResearch.Solvers.fourier_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function_5;
    ε_2 = -1,
    ε_3 = 0,
    filtration_flag = false,
    filtration_step = 2,
    filtration_factor = 2^(1/10),
)
plot(x, abs.(U))
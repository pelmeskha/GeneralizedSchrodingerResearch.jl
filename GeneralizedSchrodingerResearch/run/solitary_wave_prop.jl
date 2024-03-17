# include("run/solitary_wave_prop.jl")
using GeneralizedSchrodingerResearch

model=3

b_1 = 1
k = 0.5
omega = 0.2
theta_0 = 0
z_0 = -10
tspan = (0,1000)
xspan = (-100, 100)
h=0.25
tau=h^2

initial_function = (x) -> GeneralizedSchrodingerResearch.Initials.NSE_soliton(x, b_1, k, omega, theta_0, z_0)
x, U = GeneralizedSchrodingerResearch.Solvers.fourier_solve(
    tspan,
    xspan,
    tau,
    h,
    initial_function;
    ε_2 = -1,
    ε_3 = 0,
)
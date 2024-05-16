# include("run/validation.jl")
using Revise
using GeneralizedSchrodingerResearch.Solvers: solve
using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_5_soliton, NSE_soliton
using GeneralizedSchrodingerResearch.Utilities: construct_approximate_NSE_5_solution
using Plots, Statistics

k = 0.15
ω = 0.4

theta_0 = 0.0
z_0 = 0.0

tspan = (0.0,25.0)
xspan = (-80.0, 80.0)

h_list=[0.25, 0.2, 0.1, 0.05]
tau_list=h_list.^2

FD_error=zeros(length(h_list))
Fourier_error=zeros(length(h_list))
FD_time=zeros(length(h_list))
Fourier_time=zeros(length(h_list))

initial_function_3 = (x) -> NSE_3_soliton(x, 0.0, k, ω, theta_0, z_0)
analytical_solution_3 = (x, t) -> NSE_3_soliton(x, t, k, ω, theta_0, z_0; cycle=true, L=xspan[2]-xspan[1], c=2*k)

for i in eachindex(h_list)
    local h=h_list[i]
    local tau=tau_list[i]
    t1=time_ns()
    _, _, _, (_, _), tolerance, (_, _) = solve(
        tspan,
        xspan,
        tau,
        h,
        initial_function_3;
        method = "finite_difference",
        ε_2 = 0.0,
        ε_3 = 0.0,
        tolerance_flag = true,
        analytical_solution = analytical_solution_3,
        integrals_flag = false,
    )
    t2=time_ns()
    FD_time[i]=round((t2-t1)/1e9,digits=2)
    FD_error[i]=maximum(tolerance)
    
    t1=time_ns()
    _, _, _, (_, _), tolerance, (_, _) = solve(
        tspan,
        xspan,
        tau,
        h,
        initial_function_3;
        method = "fourier",
        ε_2 = 0.0,
        ε_3 = 0.0,
        tolerance_flag = true,
        analytical_solution = analytical_solution_3,
        integrals_flag = false,
    )
    t2=time_ns()
    Fourier_time[i]=round((t2-t1)/1e9,digits=2)
    Fourier_error[i]=maximum(tolerance)
end
xlims = (minimum(log.(h_list)),maximum(log.(h_list)))
xgap = xlims[2]-xlims[1]
xlims = (xlims[1]-xgap*0.15, xlims[2]+xgap*0.2)
ylims = (minimum(log.(Fourier_error)),maximum(log.(FD_error)))
ygap = ylims[2]-ylims[1]
ylims = (ylims[1]-ygap*0.1, ylims[2]+ygap*0.1)
plot_1 = plot(log.(h_list),log.(FD_error); xlims=xlims, ylims=ylims, label="метод конечных разностей", line=(:path,:solid,:black,2), marker=(:circle,4,:black,:black))
for i in eachindex(h_list)
    annotate!(log(h_list[i])-0.12, log(FD_error[i])+0.12, text("$(FD_time[i]) сек", :left, 10))
end
plot!(log.(h_list),log.(Fourier_error); label="метод Фурье", lw=2, line=(:path,:dash,:black,2), marker=(:d,4,:black,:gray))
for i in eachindex(h_list)
    annotate!(log(h_list[i])+0.015, log(Fourier_error[i])-0.05, text("$(Fourier_time[i]) сек", :left, 10))
end

plot!(legend=:topleft, tickfontsize=10, legendfontsize=8, yguidefontrotation=0.0)
xlabel!("log(h)")
ylabel!("логарифм относительной ошибки")
savefig(plot_1, "run/step_error.png")

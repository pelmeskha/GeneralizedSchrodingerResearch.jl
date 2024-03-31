module GeneralizedSchrodingerResearch

include("analytical_solutions/analytical_solutions.jl")
export AnalyticalSolutions

include("utilities/utilities.jl")
export Utilities

include("solvers/solvers.jl")
export Solvers

end
module GeneralizedSchrodingerResearch

include("base_utilities/base_utilities.jl")
export BaseUtilities

include("analytical_solutions/analytical_solutions.jl")
export AnalyticalSolutions

include("utilities/utilities.jl")
export Utilities

include("solvers/solvers.jl")
export Solvers

end
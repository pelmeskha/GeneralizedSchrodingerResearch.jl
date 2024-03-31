module Utilities

using GeneralizedSchrodingerResearch.AnalyticalSolutions: NSE_5_soliton

include("integrals.jl")
export integral_1, integral_2

include("filtration.jl")
export filtration

include("construct_approximate_solutions.jl")
export construct_approximate_NSE_5_solution

end
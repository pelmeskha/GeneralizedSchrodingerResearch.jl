module AnalyticalSolutions

using Interpolations
using GeneralizedSchrodingerResearch.BaseUtilities: ε2_ε3_to_M0_M1

include("NSE_3_soliton.jl")
export NSE_3_soliton

include("NSE_3_5_soliton.jl")
export NSE_3_5_soliton

include("NSE_3_5_7_soliton.jl")
export precompile_NSE_3_5_7_soliton, NSE_3_5_7_soliton

end
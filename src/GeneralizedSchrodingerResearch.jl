module GeneralizedSchrodingerResearch

include("base_utilities/base_utilities.jl")
export BaseUtilities

include("utilities/Utilities.jl")
export Utilities

include("plotting/Plotting.jl")
export Plotting

include("../Coupled_GL/src/CoupledGL.jl")
export CoupledGL

include("../GNSE_3_5_7/src/GNSE_3_5_7.jl")
export GNSE_3_5_7

end

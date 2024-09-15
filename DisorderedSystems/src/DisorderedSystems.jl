"""
Package containing the modules, functions and methods for 
the project 'Statistical Physics of Disordered Systems'.
"""
module DisorderedSystems

include("DataHandling/DataHandling.jl")
include("BeliefPropagation/BeliefPropagation.jl")
include("RandomMatrixTheory/RandomMatrixTheory.jl")
include("CavityMethod/CavityMethod.jl")
include("ProjectPlots/ProjectPlots.jl")

end

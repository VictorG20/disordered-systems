"""
Package containing the relevant functionality for the project Statistical Physics of Disordered Systems.
"""
module DisorderedSystems

include("BeliefPropagation/BeliefPropagation.jl")
export BeliefPropagation

include("RandomMatrixTheory/RandomMatrixTheory.jl")
export RandomMatrixTheory

include("ProjectPlots/ProjectPlots.jl")
export ProjectPlots

end # module DisorderedSystems

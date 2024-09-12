"""
Belief Propagation for the 1-dimensional Ising chain.
"""
module BeliefPropagation

include("cavities.jl")
include("marginals.jl")
include("magnetization.jl")
export getCavityParameters, getMeanMagnetization

end
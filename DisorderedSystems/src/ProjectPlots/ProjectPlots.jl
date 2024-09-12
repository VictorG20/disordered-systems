"""
Plotting functions for the project.
"""
module ProjectPlots

using LaTeXStrings, Plots
export display, savefig

# Belief Propagation plots. 
include("cavity_parameters.jl")
export plotCavityParameters

using DisorderedSystems.BeliefPropagation
include("magnetization.jl")
export plotMeanMagnetization

end # module ProjectPlots

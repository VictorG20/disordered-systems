"""
Plotting functions for the project.
"""
module ProjectPlots

using LaTeXStrings, Plots
export display, savefig

# Belief Propagation plots. 
include("cavity_parameters.jl")
export plotCavityParameters

using DisorderedSystems.BeliefPropagation: getMeanMagnetization
include("magnetization.jl")
export plotMeanMagnetization

include("spectral_density.jl")
export plotSpectralDensity

end # module ProjectPlots

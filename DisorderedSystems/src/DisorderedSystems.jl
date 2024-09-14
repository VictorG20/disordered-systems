"""
Package containing the relevant functionality for the project Statistical Physics of Disordered Systems.
"""
module DisorderedSystems

include("BeliefPropagation/BeliefPropagation.jl")
include("RandomMatrixTheory/RandomMatrixTheory.jl")
include("ProjectPlots/ProjectPlots.jl")

using HDF5
include("DataHandling/save.jl")
export saveData2h5

include("DataHandling/load.jl")
export loadDataFromH5

include("DataHandling/view.jl")
export viewH5FileContents

end # module DisorderedSystems

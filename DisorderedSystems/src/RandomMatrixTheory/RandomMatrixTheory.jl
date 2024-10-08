"""
This module provides the required methods to sample Gaussian Wigner Matrices
and the module to implement the cavity method.
"""
module RandomMatrixTheory

using Distributions
using LightGraphs
using LinearAlgebra
using Random

include("sample.jl")
export sampleGaussianWignerMatrix

using DisorderedSystems.DataHandling: loadDataFromH5, saveData2h5
include("direct_diagonalization.jl")
export direct_diagonalization

end
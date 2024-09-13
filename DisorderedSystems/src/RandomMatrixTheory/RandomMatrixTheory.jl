"""
This module provides the required methods to sample Gaussian Wigner Matrices
and the module to implement the cavity method.
"""
module RandomMatrixTheory

using Distributions
using LightGraphs
using LinearAlgebra
using Random
using SparseArrays

include("eigenvalues.jl")
include("sample.jl")
export sampleEigenvalues

end
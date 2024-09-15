"""
Implementation of the cavity method to estimate the 
spectral density of Gaussian Wigner Matrices.
"""
module CavityMethod
    using LinearAlgebra
    export Matrix, Symmetric, eigvals

    using SparseArrays
    export SparseMatrixCSC
    
end
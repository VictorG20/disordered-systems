
"""
    getEigenvalues(matrix::Symmetric)

Get the eigenvalues of a symmetric dense matrix.

Wrapper for the `LinearAlgebra.eigvals` method for symmetric dense matrices.
"""
function getEigenvalues(matrix::Symmetric)
    λ_values = eigvals(matrix)
    return λ_values
end


"""
    getEigenvalues(matrix::SparseMatrixCSC)

Get the eigenvalues of a symmetric sparse matrix by first turning
it into a symmetric dense matrix.
"""
function getEigenvalues(matrix::SparseMatrixCSC)
    # Make full symmetric matrix and compute eigenvalues.
    full_matrix = Symmetric(Matrix(matrix))  
    λ_values = getEigenvalues(full_matrix)
    return λ_values
end

"""
    sampleGaussianWignerMatrix(linear_size::Int, mean_connectivity::Int; 
    mean::Real = 0, std::Real = 1/sqrt(mean_connectivity), sparse::Bool = false)

Sample a real Gaussian Wigner Matrix.

The matrix is generated from the element-wise product of the adjacency matrix of 
a random regular graph of the given size and mean connectivity, with a symmetric
matrix whose entries are sampled from a normal distribution, defined by the mean
and standard deviation.

# Keyword Arguments
- `mean::Real=0`: Mean of the Gaussian distribution.
- `std::Real=1/\\sqrt{mean_connectivity}`: Standard deviation of the Gaussian distribution.
- `sparse::Bool=false`: If true, return the matrix as type `SparseMatrixCSC`.
"""
function sampleGaussianWignerMatrix(linear_size::Int, mean_connectivity::Int; 
    mean::Real = 0, std::Real = 1/sqrt(mean_connectivity), sparse::Bool = false)

    rrg = random_regular_graph(linear_size, mean_connectivity)
    A = adjacency_matrix(rrg)
    weights_distro = Normal(mean, std)
    J = Symmetric(rand(weights_distro, (linear_size, linear_size)))
    M = A .* J

    # If not sparse, return as a full, symmetric, matrix.
    M = sparse ? (M .* A) : Symmetric(M)

    return M
end
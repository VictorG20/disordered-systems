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


"""
    sampleEigenvalues(linear_size::Int, mean_connectivity::Int, num_samples::Int; 
    mean::Real = 0, std::Real = 1/sqrt(mean_connectivity))

Sample eigenvalues of real Gaussian Wigner matrices.

`num_samples` real Gaussian Wigner matrices of `linear_size` and `mean_connectivity`
are generated and their eigenvalues computed via direct diagonalization. The output 
is a matrix of eigenvalues where each column corresponds to a single sample.

# Keyword Arguments
- `mean::Real=0`: Mean of the Gaussian distribution.
- `std::Real=1/\\sqrt{mean_connectivity}`: Standard deviation of the Gaussian distribution.
"""
function sampleEigenvalues(linear_size::Int, mean_connectivity::Int, num_samples::Int; 
    mean::Real = 0, std::Real = 1/sqrt(mean_connectivity))

    λ_samples = zeros(Float64, linear_size, num_samples)

    for sample in 1:num_samples
        M = sampleGaussianWignerMatrix(linear_size, mean_connectivity; mean=mean, std=std)
        λ_samples[:, sample] = eigvals(M)
    end

    return λ_samples
end
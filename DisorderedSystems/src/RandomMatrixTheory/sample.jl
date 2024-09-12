"""
    sampleGaussianWignerMatrix(N::Int, c::Int; μ::Real = 0, σ::Real = 1/sqrt(c))

Sample a (real) Gaussian Wigner Matrix of size N x N, mean connectivity c.
The Gaussian distribution has by default a mean μ = 0 and a standard deviation
of σ = 1/sqrt(c).
"""
function sampleGaussianWignerMatrix(N::Int, c::Int; μ::Real = 0, σ::Real = 1/sqrt(c))

    rrg = random_regular_graph(N, c)
    A = adjacency_matrix(rrg)
    weights_distro = Normal(μ, σ)
    J = Symmetric(rand(weights_distro, (N, N)))
    M = A .* J
    M = M .* A  # This passes the matrix as a sparse matrix.
    return M
end
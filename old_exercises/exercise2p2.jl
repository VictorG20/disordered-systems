using SparseArrays
using Random

using DisorderedSystems.CavityMethod
using DisorderedSystems.RandomMatrixTheory: sampleGaussianWignerMatrix


function getEigenvalues(M::SparseMatrixCSC)
    dense_matrix = Symmetric(Matrix(M))
    eigenvalues = eigvals(dense_matrix)
    return eigenvalues
end


function getCavitiesMatrix(M::SparseMatrixCSC{T}) where {T}

    N = size(M)[1]
    c = length(M[:, 1].nzind)

    rows = []
    cols = []
    vals = T[]
    
    for j in 1:N
        j_neighbors = M[:, j].nzind
        for (k_idx, k) in enumerate(j_neighbors)
            row_idx = (j - 1) * c + k_idx
            k_neighbors = M[:, k].nzind
            for (l_idx, l) in enumerate(k_neighbors)
                (l == j) && continue
                col_idx = (k - 1) * c + l_idx  # Index of cavity precision in vector.
                push!(rows, row_idx)
                push!(cols, col_idx)
                push!(vals, M[k, l]^2)
            end
        end
    end

    M_cavities = sparse(rows, cols, vals, N*c, N*c)

    return M_cavities
end


function getCavityPrecisions(M, z)
    iterations = 1000
    check_every = 100
    abs_tol = 1e-6

    N = size(M)[1]
    c = length(M[:, 1].nzind)

    M_cavities = getCavitiesMatrix(M)
    ω_cavities = rand(Complex{Float64}, N*c)

    tol_counter = 0

    for iteration in 1:iterations
        
        tol_counter += 1

        if tol_counter != check_every
            ω_cavities = z .+ M_cavities * (1. ./ ω_cavities)
            continue
        end

        tol_counter = 0
        new_omegas = z .+ M_cavities * (1. ./ ω_cavities)
        Δω = new_omegas - ω_cavities
        maximum_Δω = max(maximum(abs.(real(Δω))), maximum(abs.(imag(Δω))))
        ω_cavities = new_omegas

        if maximum_Δω < abs_tol
            println("Reached convergence in $(iteration) iterations")
            break
        end
    end
    
    ω_marginal = z .+ vec(sum(reshape((M.nzval .^2) ./ ω_cavities, (N, c)), dims=2))
    ρ = sum(imag.(im ./ ω_marginal)) / (π * N)
end


function main()
    N = 2^11
    c = 3
    ϵ = 0.005
    λ_range = range(0, 3, step=0.01)
    samples = 1

    Random.seed!(42)
    M = sampleGaussianWignerMatrix(N, c; sparse=true)
    
    
    # λ_values = getEigenvalues(M)
    println(λ_range[50])
    z_range = (im * λ_range) .+ ϵ

    getCavityPrecisions(M, z_range[50])

end

main()

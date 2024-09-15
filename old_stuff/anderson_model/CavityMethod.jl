using Distributions
using LightGraphs
using Random


function get_cavity_precisions()
    # How many cavity precisions do I have
    
end


function compute_marginals(λ:: Float64, z:: Vector{ComplexF64}, cavity_ω:: Vector{ComplexF64})
    N = length(z)
    c = length(cavity_ω) ÷ N

    # Reshape the cavities into a 2D array such that the columns 
    # correspond to cavities (j) and the rows are the neighbors k. 
    cavities_matrix = reshape(cavity_ω, N, c)

    # Sum along the rows and add 'z' to the resulting vector.
    marginals_ω = vec(sum(cavities_matrix, dims=2))
    marginals_ω = (z + marginals_ω) .+ im*λ

    return marginals_ω
end


function main(; N:: Int, c:: Int, W:: Float64, ϵ:: Float64)
    Random.seed!(42)
    λ_values = [0.]

    # Create the adjacency matrix.
    RRG = random_regular_digraph(N, c)
    A = adjacency_matrix(RRG)

    # Sample the on-site energies.
    energy_distro = Uniform(-W/2, W/2)
    energies = rand(energy_distro, N)

    z = ϵ .- im*energies  # Constant term in the equations for the cavity and marginal precisions.

    # Initialize array to store spectral_density.
    ρ_values = Array{Float64}(undef, length(λ_values))

    # Initialize cavity precisions.
    cavity_ω = rand(ComplexF64, N*c)

    for (k, λ) in enumerate(λ_values)
        # Lead cavity precisions to fixed point with given accuracy.

        # Get marginals.
        marginal_ω = compute_marginals(λ, z, cavity_ω)

        # Compute spectral density and store value.
        ρ_values[k] = sum(imag.(im ./ marginal_ω)) / (π * N)
    end
    
    return ρ_values
end

main(N=8, c=3, W=0.3, ϵ=1e-3)


# TODO: Check creating ids array vs. simply accessing the neighbors of the adjacency matrix.
# TODO: Check if it's better (faster with same accuracy) to track changes of real and imaginary parts separately.

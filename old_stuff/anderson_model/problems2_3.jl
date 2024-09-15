include("PopulationDynamics.jl")

using Distributions
using LinearAlgebra
using LightGraphs
using Random


function direct_diagonalization(E, A)
    H = Symmetric(Matrix(Diagonal(E) - A))
    λ_values = eigvals(H)
    
    return λ_values
end


function main()
    Random.seed!(42)  # Fix seed for reproducibility.

    # System parameters
    N:: Int = 2^10  # System size.
    c:: Int = 3  # Mean connectivity.
    W:: Float64 = 0.3  # Disorder.

    ###########################
    # Direct diagonalization. #
    ###########################
    RRG = random_regular_graph(N, c)
    A = adjacency_matrix(RRG)
    E = rand(Uniform(-W/2, W/2), N)
    
    println("Starting direct diagonalization")
    λ_values = direct_diagonalization(E, A)

    ########################
    # Population Dynamics. #
    ########################

    println("Starting Population Dynamics")
    # Parameters.
    sample_size::Int = 100  # Number of marginal precisions for measurement.
    ϵ::Float64 = 1e-3  # Imaginary part to eigenvalues.
    Np:: Int = 10^3  # Population size / number of cavity precisions.
    max_sweeps:: Int = 100  # Maximum number of PD sweeps to reach equilibrium.
    
    λ_range = range(0, 3, length=150)  # Range of λ_values to probe.

    PDSyste = PDSystem(sample_size, Np, c, ϵ, W, max_sweeps)
    pd_density = spectraldensity(PDSyste, λ_range)

    return pd_density
end

main()

# TODO: Read parameters from file (.json, .dat or something like that) and save results in .h5 or similar.
# TODO: Create a single plot with histogram, cavity method and population dynamics results.
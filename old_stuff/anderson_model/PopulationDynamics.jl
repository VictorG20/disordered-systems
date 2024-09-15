# module PopulationDynamics

using Distributions
using ProgressMeter
using Random


struct PDSystem
    sample_size:: Int
    population_size:: Int  # Population size
    connectivity:: Int  # Mean connectivity
    epsilon:: Float64  # Imaginary part to eigenvalues.
    disorder:: Float64  # Disorder for on-site energies.
    equilibrium_sweeps:: Int  # Maximum number of sweeps for equilibrium.
end


"""
    pdsweep!(Np:: Int, c:: Int, ϵ:: Float64, energy_distro, cavities:: Vector{ComplexF64}, λ:: Float64)

TBW
"""
function pdsweep!(Np:: Int, c:: Int, ϵ:: Float64, energy_distro, cavities:: Vector{ComplexF64}, λ:: Float64)
    # Random order to change every element of the population
    rand_ids = sample(1:Np, Np; replace=false)

    for idx in rand_ids
        # Sample c-1 cavity precisions, energy, and update element of population.
        ω_sample = sample(cavities, c-1, replace=false)
        E = rand(energy_distro)

        cavities[idx] = im*(λ - im*ϵ - E) + sum(1 ./ ω_sample)
    end
end


"""
    samplemarginal(c:: Int, ϵ:: Float64, energy_distro, cavities:: Vector{ComplexF64}, λ:: Float64)

TBW
"""
function samplemarginal(c:: Int, ϵ:: Float64, energy_distro, cavities:: Vector{ComplexF64}, λ:: Float64)
    # Sample c cavity precisions, energy, and compute marginal precision.
    cavity_sample = sample(cavities, c, replace=false)
    E = rand(energy_distro)

    marginal_ω = im*(λ - im*ϵ - E) + sum(1 ./ cavity_sample)

    return marginal_ω
end


"""
    spectraldensity(system:: PDSystem, λ_values; fix_seed:: Bool = false)

TBW
"""
function spectraldensity(system:: PDSystem, λ_values; fix_seed:: Bool = false)
    fix_seed && Random.seed!(42)  # Fix the seed for reproducibility.

    N = system.sample_size
    Np = system.population_size
    W = system.disorder
    c = system.connectivity
    ϵ = system.epsilon
    equilibrium_sweeps = system.equilibrium_sweeps

    energy_distro = Uniform(-W/2, W/2)  
    sweep!(cavities:: Vector{ComplexF64}, λ::Float64) = pdsweep!(Np, c, ϵ, energy_distro, cavities, λ)
    getmarginal(cavities:: Vector{ComplexF64}, λ::Float64) = samplemarginal(c, ϵ, energy_distro, cavities, λ)

    # Initialize complex population of cavity precisions and array for spectral density values.
    cavities_ω = rand(Complex{Float64}, Np)
    ρ_values = Array{Float64}(undef, length(λ_values))

    @showprogress for (k, λ) in enumerate(λ_values)

        # Equilibrate cavities.
        for _ in 1:equilibrium_sweeps
            sweep!(cavities_ω, λ)
        end

        # Sample N marginal precisions, performing a pd_sweep after each sample.
        marginals_ω = zeros(ComplexF64, N)
        for j in 1:N
            marginals_ω[j] = getmarginal(cavities_ω, λ)
            sweep!(cavities_ω, λ)
        end
        
        ρ_values[k] = sum(imag(im ./ marginals_ω)) / (π * N)
    end

    return ρ_values
end

# end;  # End of module.

# TODO: Document functions.
# TODO: Devise method to check equilibrium of the cavity precisions.
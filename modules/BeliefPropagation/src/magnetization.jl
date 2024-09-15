"""
    getMeanMagnetization(marginal_probs::Vector{Float64}, N::Int; up::Bool=true)

Compute the mean magnetization of a chain of 'N' nodes via its marginal probabilities.
If 'up', the marginal probabilities correspond to P(σ = +1). If less than 'N' marginal
probabilities are given, it is assumed that the chain is symmetric at its midpoint and
the remaining marginal probabilities are computed accordingly.
"""
function getMeanMagnetization(marginal_probs::Vector{Float64}, N::Int; up::Bool=true)    
    
    S = length(marginal_probs)

    if S > N
        error("A chain with $(N) nodes has at most $(N) 
            distinct marginals but $(S) were given.")
    elseif S == 0
        error("Marginal probabilities vector is empty.")
    end

    # Compute avg. σ = P(σ=+1) - P(σ=-1).
    avgs_σ = 2*marginal_probs .- 1 
    (!up) && (avgs_σ *= -1)  # If marginals are P(σ = -1).

    if S == N
        # We have all the averages.
        mean_mag = sum(avgs_σ) / N
    elseif S < N/2
        # Fixed point was achieved. Complete missing values.
        mean_mag = (2*sum(avgs_σ) + (N - 2*S)*avgs_σ[end]) / N
    else
        mean_mag = (sum(avgs_σ) + sum(avgs_σ[1:N-S])) / N
    end
    return mean_mag
end


"""
    getMeanMagnetization(β::Real, h::Real)

Get the mean magnetization at an inverse thermal energy 'β' and external
magnetic field 'h' for the Ising chain in the thermodynamic limit. This 
corresponds to the analytical expression for the mean magnetization.
"""
function getMeanMagnetization(β::Real, h::Real)
    x = exp(2. * β) * sinh(β * h)
    M = x / sqrt(1. + x^2)
    return M
end


"""
    getMeanMagnetization(N::Int, β_values::Vector{Float64}, 
    h::Real, init_prob::Real, up::Bool)

Get the mean magnetization for a series of inverse thermal energies β and external magnetic
field h, for an Ising chain of N nodes. This method computes first the cavity parameters 
of the chain, for which an initial probability, corresponding to the cavity probability at
the one of the ends of the chain, is required. The parameter 'up' refers to the orientation
of the spin for the probability given: if true, the initial probability is for P(σ = +1).
"""
function getMeanMagnetization(N::Int, β_values::Vector{Float64}, 
    h::Real, init_prob::Real, up::Bool)

    m_values = Array{Float64}(undef, length(β_values))

    for (k, β) in enumerate(β_values)
        cavity_params = getCavityParameters(N, β, h, init_prob, up=up)
        cavity_probs = getCavityProbability.(cavity_params, β, up=true)
        margin_probs = getMarginalProbabilities(cavity_probs, N, β, h)
        m_values[k] = getMeanMagnetization(margin_probs, N)
    end

    return m_values
end
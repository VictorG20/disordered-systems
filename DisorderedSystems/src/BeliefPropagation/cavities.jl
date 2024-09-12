
"""
    getCavityRecursively(ω::Float64, β::Real, h::Real)

Compute the next cavity parameter in a chain via the recursive relation.
This relation holds for both left and right messages.
"""
function getCavityRecursively(ω::Float64, β::Real, h::Real)
    return h + (1/(2 * β)) * log(cosh(β*(ω + 1)) / cosh(β*(ω - 1)))
end


"""
    getCavityFromProbability(probability::Real, β::Real; up=true)

Compute the cavity parameter from the cavity probability. 
If 'up' is true, 'probability' is P(σ= +1), else P(σ=-1).
"""
function getCavityFromProbability(probability::Real, β::Real; up=true)

    ω = log(probability / (1 - probability)) / (2 * β)
    (!up) && (ω *= -1)  # The probability is for P(σ=-1), change sign.

    return Float64(ω)
end


"""
    getCavityParameters(N::Int, β::Real, h::Real, init_prob::Real; 
    up::Bool=true, abstol::Float64=1e-10, fill::Bool=false, print_fixed::Bool=false)

Compute the 'N-1' cavity parameters of either the left or right messages of a chain of 
'N' spins with inverse thermal energy 'β' and external magnetic field strength 'h'. The 
first cavity parameter is computed from its cavity probability 'init_prob' corresponding 
to P(σ = +1) if 'up' is true (default), or P(σ = -1) if false; the remaining parameters 
are computed recursively. The recursive computation stops if a fixed point is achieved 
with 'abstol'. A vector with the cavity parameters up to the fixed point is returned or,
if 'fill', an 'N-1' vector is returned with the fixed point at every subsequent place.
"""
function getCavityParameters(N::Int, β::Real, h::Real, init_prob::Real; 
    up::Bool=true, abstol::Float64=1e-10, fill::Bool=false, print_fixed::Bool=false)

    # A chain with N nodes has N-1 cavity parameters.
    cavity_params = Array{Float64}(undef, N-1)
    cavity_params[1] = getCavityFromProbability(init_prob, β, up=up)

    # Compute the remaining cavity parameters recursively until fixed point.
    idx = 0  # Save last index of array in case fixed point is reached.
    for i in 2:N-1
        # Get next cavity parameter and compute difference.
        cavity_params[i] = getCavityRecursively(cavity_params[i-1], β, h)
        Δω = abs(cavity_params[i] - cavity_params[i-1])

        # If the difference is below the tolerance, stop.
        if Δω < abstol
            idx = i-1  # Save index of last cavity parameter required.
            print_fixed && println("    Fixed point value achieved: ", cavity_params[idx], 
                                   " with abstol: $(abstol) after $(idx) iterations.")
            break
        end
    end

    # If fixed point was reached...
    if fill && (0 < idx < N - 1)
        # Fill remaining cavity parameters with fixed point.
        cavity_params[idx+1:end] .= cavity_params[idx]
    elseif (0 < idx < N - 1)
        # Drop undefined values.
        cavity_params = cavity_params[1:idx]
    end

    return cavity_params
end


"""
    getCavityProbability(ω::Float64, β::Real; up::Bool=true)

Compute the cavity probability from the cavity parameter and an inverse
thermal energy 'β'. If 'up', the returned probability is for P(σ = +1).
"""
function getCavityProbability(ω::Float64, β::Real; up::Bool=true)
    σ = up ? 1. : -1.
    cavity_prob = (1. + σ * tanh(β * ω)) / 2.
    return cavity_prob
end
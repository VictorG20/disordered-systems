"""
    getMarginalProbabilities(cavity_probs::Vector{Float64}, N::Int, β::Real, h::Real; 
    up_in::Bool=true, up_out::Bool=true)
    
Compute the marginal probabilities of a chain of 'N' spins, at inverse thermal energy 'β'
and external magnetic field 'h', from the (distinct) cavity probabilities. The parameters
'up_in' and 'up_out' specify if the input/output probabilities are for P(σ = +1).

It is assumed that, if the length of 'cavity_probs' is less than 'N', the fixed point was
achieved. From this, only the distinct marginals necessary to know all of them are obtained.
"""
function getMarginalProbabilities(cavity_probs::Vector{Float64}, N::Int, β::Real, h::Real; 
    up_in::Bool=true, up_out::Bool=true)
    
    S = length(cavity_probs)

    if S > N - 1
        error("A chain of $(N) nodes has at most $(N-1) distinct 
            cavity parameters but $(S) were given.")
    end

    marginal_prob_up = 0.  # P(σi = +1).
    marginal_prob_dn = 0.  # P(σi = -1).

    for σi in [1., -1.]
        if up_in
            left_message = (exp(β*σi) * cavity_probs) + (exp(-β*σi) * (1. .- cavity_probs))
        else
            left_message = (exp(β*σi) * (1. .- cavity_probs)) + (exp(-β*σi) * cavity_probs)
        end

        marginal_prob = zeros(Float64, S + 1)
        marginal_prob[1] = left_message[end]  # End of left message == end of right message.

        if S < (N-1)/2
            # Relevant part of right message is at fixed point.
            marginal_prob[2:end] = left_message * left_message[end]
        elseif S < N - 1
            # Reached fixed point. Relevant right message is mixed.
            # There are N-S values with the fixed point, one goes 
            # to the first place of the marginals, N-S-1 remain.
            marginal_prob[2:N-S] = left_message[1:N-S-1] * left_message[end]
            marginal_prob[N-S+1:S+1] = left_message[N-S:S] .* reverse(left_message[N-S-1:S-1])
        else
            # Fixed point wasn't reached. Compute all marginals.
            marginal_prob[2:N-1] = left_message[1:N-2] .* (reverse(left_message))[2:N-1]
            marginal_prob[N] = left_message[N-1]
        end

        marginal_prob *= exp(β*h*σi)
        
        if σi == 1
            marginal_prob_up = marginal_prob
        else
            marginal_prob_dn = marginal_prob
        end
    end

    # Normalize the probabilities.
    norm_factor = marginal_prob_up + marginal_prob_dn

    (up_out) && return marginal_prob_up ./ norm_factor

    return marginal_prob_dn ./ norm_factor
end

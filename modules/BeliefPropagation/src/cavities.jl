@doc raw"""
    get_cavity_parameter(ω::T, β::T, h::T) where T <: Real

Get recursively the next cavity parameter in the chain.

The recursive relation between the cavity parameters is given by
```math
\omega_{i-1}^{(i)} = h + \frac{1}{2\beta} \log{\left( 
    \frac{
    \cosh{[\beta(\omega_{i-2}^{(i-1)} + 1)]}
    }{
    \cosh{[\beta(\omega_{i-2}^{(i-1)} - 1)]}
    }
    \right)}
```
"""
function get_cavity_parameter(ω::T, β::T, h::T) where T <: Real
    likelihood = cosh(β * (ω + 1)) / cosh(β * (ω - 1))
    next_omega = h + (1 / (2 * β)) * log(likelihood)
    return next_omega
end



# """
#     getCavityRecursively(ω::Float64, β::Real, h::Real)

# Compute the next cavity parameter in a chain via the recursive relation.
# This relation holds for both left and right messages.
# """
# function getCavityRecursively(ω::Float64, β::Real, h::Real)
#     return h + (1/(2 * β)) * log(cosh(β*(ω + 1)) / cosh(β*(ω - 1)))
# end


# """
#     getCavityFromProbability(probability::Real, β::Real; up=true)

# Compute the cavity parameter from the cavity probability. 
# If 'up' is true, 'probability' is P(σ= +1), else P(σ=-1).
# """
# function getCavityFromProbability(probability::Real, β::Real; up=true)

#     ω = log(probability / (1 - probability)) / (2 * β)
#     (!up) && (ω *= -1)  # The probability is for P(σ=-1), change sign.

#     return Float64(ω)
# end
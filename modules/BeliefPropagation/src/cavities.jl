"""
    getCavityRecursively(ω::Real, β::Real, h::Real)

Get the next cavity parameter for a chain at inverse temperature β and 
external magnetic field h, using the recursive relation between them.
"""
function getCavityRecursively(ω::Real, β::Real, h::Real)
    likelihood = cosh(β * (ω + 1)) / cosh(β * (ω - 1))
    next_omega = h + (1 / (2 * β)) * log(likelihood)
    return next_omega
end


"""
    getCavityFromProbability(probability::Real, β::Real; up=true)

Compute the cavity parameter from the cavity probability. 
If 'up' is true, 'probability' is P(σ= +1), else P(σ=-1).
"""
function getCavityFromProbability(probability::Real, β::Real; up=true)
    likelihood = probability / (1 - probability)
    ω = log(likelihood) / (2 * β)
    (!up) && (ω *= -1)  # The probability is for P(σ=-1), change sign.
    return ω
end
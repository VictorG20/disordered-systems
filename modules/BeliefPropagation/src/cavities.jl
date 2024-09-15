"""
    get_cavity_parameter(ω::T, β::T, h::T) where T <: Real

Get recursively the next cavity parameter for a chain at inverse temperature β and
external magnetic field h.
"""
function get_cavity_parameter(ω::T, β::T, h::T) where T <: Real
    likelihood = cosh(β * (ω + 1)) / cosh(β * (ω - 1))
    next_omega = h + (1 / (2 * β)) * log(likelihood)
    return next_omega
end

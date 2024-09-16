"""
Belief propagation algorithm methods for solving 
the one-dimensional Ising model (Ising chain).
"""
module BeliefPropagation


@kwdef struct Chain
    N::Int   # Number of nodes in the chain.
    h::Real  # External magnetic field.
    J::Real  # Coupling constant
end

# Base.show(io::IO, ::MIME"text/plain", chain::Chain) = print(
#     io, "Examplary instance of MyType\n", chain.N, " ± ", chain.h)
Base.show(io::IO, chain::Chain) = print(
    io, "Chain(N=$(chain.N), h=$(chain.h), J=$(chain.J))")

export Chain


include("cavities.jl")
export get_cavity_parameter

end # module BeliefPropagation

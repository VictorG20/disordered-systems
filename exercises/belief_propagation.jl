using BeliefPropagation


function cavity_vs_position()
    chain1 = Chain(N=10, h=1//2, J=1.)
    chain2 = Chain(N=20, h=0.7, J=1)
    println(chain1)
    println(typeof([chain1, chain2]))
end


cavity_vs_position()
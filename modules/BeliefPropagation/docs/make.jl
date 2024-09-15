using Documenter
using BeliefPropagation

makedocs(sitename="BeliefPropagation")

deploydocs(; repo = "github.com/VictorG20/disordered-systems.git",
             dirname = "BeliefPropagation",
             versions = nothing)
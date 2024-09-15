using Documenter
using BeliefPropagation

makedocs(sitename="BeliefPropagation")

deploydocs(; repo = "github.com/VictorG20/disordered-systems.git",
             dirname="BeliefPropagation",
             tag_prefix="modules/BeliefPropagation-",
             # ...any additional kwargs
             )
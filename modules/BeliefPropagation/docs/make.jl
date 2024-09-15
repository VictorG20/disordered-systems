using Documenter
using BeliefPropagation

makedocs(sitename="BeliefPropagation")

deploydocs(; repo = "github.com/VictorG20/disordered-systems.git",
             dirname = "BeliefPropagation",
             versions = nothing,
            #  devbranch = "main",
            #  deploy_config = Documenter.GitHubActions("github.com/VictorG20/disordered-systems.git", "workflow_dispatch", "main"),
            #  tag_prefix = "BeliefPropagation-",
        )
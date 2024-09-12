"""
    plotCavityParameters(cavities_list::Vector{Vector{Float64}}, probabilities::Vector{Float64};my_markerstrokewidth::Float64=2.7)

Plot cavity parameters vs. cavity positions.

The vector 'cavities_list' should contain, in each entry, the cavity parameters of the chain.
If the chain is of length 'N', there are 'N-1' cavity parameters.
"""
function plotCavityParameters(cavities_list::Vector{Vector{Float64}}, probabilities::Vector{Float64}; my_markerstrokewidth::Float64=2.7)
    
    # ----------------------------------------------------------------------
    # General plot attributes. 

    kw = (;
    xlabel=L"Cavity position $i$",
    ylabel="\$ \\omega^{(i)}_{i-1} \$",
    xminorgrid=true,
    yminorgrid=true,
    xminorticks=5,
    yminorticks=5,
    xlabelfontsize=12,
    ylabelfontsize=14,  # Default value is 11.
    legendfontsize=13,
    legendtitlefontsize=15,
    legend=:topright,
    legend_columns=-1,
    tickfontsize=13,
    framestyle=:box,
    dpi=400,
    size=(640, 300),
    bottommargin=4Plots.mm,
    leftmargin=3Plots.mm
    )

    N = length(cavities_list[1]) + 1  # There are N-1 cavity parameters.
    
    # Ticks in LaTeX for uniform font.
    xticks_range = (N > 15) ? (5:5:N) : (0:2:N)
    xticks_labs = ["\$ $(k) \$" for k ∈ xticks_range]
    yticks_range = -1:2
    yticks_labs = ["\$ $(k) \$" for k ∈ yticks_range]

    # ----------------------------------------------------------------------
    # Series attributes.

    marker_colors = [1, 2, :blue, :orange]
    marker_shapes = [:circle, :xcross, :vline, :hline]
    plot_labels = ["\$ $(probability) \$" for probability in probabilities]
    
    # ----------------------------------------------------------------------
    # Initiate plot with first probability.

    cavity_pos = 2:N  # Cavity positions.
    cavity_plot = scatter(
        cavity_pos, cavities_list[1], label=plot_labels[1], 
        markershape=marker_shapes[1], markerstrokewidth=my_markerstrokewidth, 
        markerstrokecolor=marker_colors[1], markersize=4
        )
    plot!(;kw...)
    plot!(xticks=(xticks_range, xticks_labs),
        yticks=(yticks_range, yticks_labs),)
    # Legend title looks fugly.
    # plot!(legendtitle="\$P^{(2)}(\\sigma_{1}=1)\$")
    # ----------------------------------------------------------------------

    (length(probabilities) == 1) && return cavity_plot
    
    # ----------------------------------------------------------------------
    # Plot the rest.
    for (k, cavity_params) in enumerate(cavities_list)
        (k==1) && continue
        scatter!(cavity_pos, cavity_params, label=plot_labels[k], 
            markershape=marker_shapes[k], 
            markercolor=marker_colors[k],
            markerstrokewidth=my_markerstrokewidth, 
            )
    end
    # ----------------------------------------------------------------------
    return cavity_plot
end
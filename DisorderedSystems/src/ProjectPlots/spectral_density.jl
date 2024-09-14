function _getBinEdges(data, bin_step, min_λ)
    min_idx = floor(Int, (minimum(data) - min_λ) / bin_step)
    max_idx = ceil(Int, (maximum(data) - min_λ) / bin_step)
    min_edge = min_λ + min_idx * bin_step
    bin_edges = range(min_edge, step=bin_step, length=max_idx - min_idx + 1)
    return bin_edges
end


function _formatTicks(histo_plot)

    # x-axis ticks.
    max_x_value = maximum(abs.(histo_plot.series_list[1].plotattributes[:x]))
    xticks_max = round(Int, max_x_value)
    xticks_range = -xticks_max:xticks_max
    xticks_label = ["\$ $(i) \$" for i ∈ xticks_range]
    plot!(xticks=(xticks_range, xticks_label))

    # y-axis ticks
    bar_heights = [y for y in histo_plot.series_list[1].plotattributes[:y] if !isnan(y)]
    y_max = maximum(bar_heights)
    yticks_range = [y for y in 0:0.1:1 if y <= 1.05*y_max]
    
    # Avoid saturating with many ticks.

    max_yticks = 6

    if length(yticks_range) > max_yticks
        yticks_max = maximum(yticks_range)
        yticks_range = range(0, yticks_max, length=max_yticks)
    end
    
    plot!(yminorticks=2)
    
    yticks_label = ["\$ $(round(y, digits=2)) \$" for y in yticks_range]
    plot!(yticks=(yticks_range, yticks_label))

    return nothing
end


function plotSpectralDensity(histos; kwargs...)
    
    if :shared_bins in keys(kwargs) && kwargs[:shared_bins]

        !(:number_bins in keys(kwargs)) && 
            error("`shared_bins` option requires giving a `number_bins` argument.")

        min_λ = minimum([minimum(histo.data) for histo in histos])
        max_λ = maximum([maximum(histo.data) for histo in histos])
        bin_step = (max_λ - floor(min_λ, digits=4)) / kwargs[:number_bins]

        # Create first plot
        histo_plot = stephist(
            histos[1].data; 
            histos[1].h_kw..., bins=_getBinEdges(histos[1].data, bin_step, min_λ)
            )

        # Create remaining plots.
        for histo in histos[2:end]
            bins=_getBinEdges(histo.data, bin_step, min_λ)
            stephist!(histo.data; histo.h_kw..., bins=bins)
        end

    else
        histo_plot = stephist(histos[1].data; histos[1].h_kw...)
        [stephist!(histo.data; histo.h_kw...) for histo in histos[2:end]]
    end

    # Wigner's semicircle law line plot.
    if :plot_wigner in keys(kwargs) && kwargs[:plot_wigner]
        λ_range = range(-2., 2., length=500)
        plot!(
            λ_range, sqrt.(4 .- λ_range.^2)/(2*π), 
            lc=:red, lw=2.5, label="\$ \\mathrm{Wigner's} \\;\\; \\mathrm{law} \$", 
            z_order=1
            )
    end

    # Ticks formatting
    (:format_tick in keys(kwargs)) && (kwargs[:format_tick]) && _formatTicks(histo_plot)

    # Legend title.
    (:legendtitle in keys(kwargs)) && plot!(; legendtitle = kwargs[:legendtitle])

    # Kwargs for the overall plot.
    (:plots_kwarg in keys(kwargs)) && plot!(;kwargs[:plots_kwarg]...)

    return histo_plot
end
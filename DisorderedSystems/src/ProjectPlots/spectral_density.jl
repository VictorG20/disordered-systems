function plotSpectralDensity(λ_series; 
    bins::Int=50, histo_lw::Union{Float64, Int} = 1.7, colors=nothing, labels=nothing,
    plot_wigner::Bool = false, kw=nothing)

    minima = [minimum(vec(λ_samples)) for λ_samples in λ_series]
    maxima = [maximum(vec(λ_samples)) for λ_samples in λ_series]
    min_λ  = minimum(minima)
    max_λ = maximum(maxima)

    bin_step = (max_λ - floor(min_λ, digits=4)) / bins

    if isnothing(labels)
        labels = [nothing for _ in 1:(size(λ_series)[1] + 1)]
    end

    if isnothing(colors)
        # Add one extra color for the plot of Wigner's semicircle law.
        colors = [k for k in 1:(size(λ_series)[1] + 1)]
    end
    
    histo = nothing
    
    for (k, (λ_samples, label, color)) in enumerate(zip(λ_series, labels, colors))
        min_idx = floor(Int, (minima[k] - min_λ) / bin_step)
        max_idx = ceil(Int, (maxima[k] - min_λ) / bin_step)
        min_edge = min_λ + min_idx * bin_step
        bin_edges = range(min_edge, step=bin_step, length=max_idx - min_idx + 1)

        if k == 1
            histo = stephist(
                vec(λ_samples), normalize=:pdf, bins=bin_edges, lw=histo_lw, 
                color=color, label=label, legendtitle=labels[end])
        else
         stephist!(
            vec(λ_samples), normalize=:pdf, bins=bin_edges, lw=histo_lw, 
            color=color, label=label)
        end
    end

    if plot_wigner
        λ_range = range(-2., 2., length=500)
        plot!(
            λ_range, sqrt.(4 .- λ_range.^2)/(2*π), 
            lc=colors[end], lw=2.5, label="\$ \\mathrm{Wigner's} \\;\\; \\mathrm{law} \$", z_order=1)
    end

    # x-axis ticks.
    xticks_max = round(Int, max(abs(min_λ), max_λ))
    xticks_range = -xticks_max:xticks_max
    xticks_label = ["\$ $(i) \$" for i ∈ xticks_range]
    plot!(xticks=(xticks_range, xticks_label))

    # y-axis ticks
    bar_heights = [y for y in histo.series_list[1].plotattributes[:y] if !isnan(y)]
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

    if !isnothing(kw)
        plot!(;kw...)
    end

    return histo
end
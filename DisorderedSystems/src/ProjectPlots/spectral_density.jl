function plotSpectralDensity(λ_series; plot_wigner::Bool = false)


    bins = 50
    bin_edges = range(-3, 3, length=bins+1)

    histo = scatterhist(vec(λ_series[1]), normalize=:pdf, bins=bin_edges, color=:red)
    scatterhist!(vec(λ_series[2]), normalize=:pdf, bins=bin_edges, color=:orange)

    if plot_wigner
        λ_range = range(-2., 2., length=500)
        plot!(
            λ_range, sqrt.(4 .- λ_range.^2)/(2*π), 
            lc=:blue, lw=2.4, label="Wigner's law", z_order=1)
    end

    display(histo)

    return
end
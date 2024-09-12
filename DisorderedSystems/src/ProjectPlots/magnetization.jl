
"""
    plotMeanMagnetization(m_values::Vector{Float64}, β_values::Vector{Float64}, 
    h::Real, N::Int)

Plot the mean magnetization values at different inverse thermal energy values
but fixed external magnetic field strength 'h'. 'N' is the length of the chain.
"""
function plotMeanMagnetization(m_values::Vector{Float64}, β_values::Vector{Float64}, 
    h::Real, N::Int)

    # Create series of β values for the exact solution.
    β_min = minimum(β_values)
    β_max = maximum(β_values)
    β_exact = exp10.(range(log10(β_min), log10(β_max), length=500))
    
    m_exact = getMeanMagnetization.(β_exact, h)


    m_plot = plot(β_exact, m_exact, label="Analytical", lw=2.5, lc=2)
    scatter!(
        β_values, m_values, 
        label="Numerical", 
        markercolor=1, 
        markersize=5, 
        markerstrokewidth=0.5, 
        markeralpha=0.85)
    
        # General plot attributes.

        h_label = trunc(Int, h)

    kw = (;
        xlabel="\$ \\beta \$",
        ylabel="\$ \\langle M \\, \\rangle \$",
        legendtitle="\$ h = $(h_label) \\quad N = $(N) \$",
        xminorgrid=true,
        yminorgrid=true,
    #     xminorticks=10,
        yminorticks=5,
        xlabelfontsize=14,
        ylabelfontsize=14,  # Default value is 11.
        legendfontsize=12,
        legendtitlefontsize=14,
        legend=:right,
    #     legend_columns=-1,
        tickfontsize=13,
        framestyle=:box,
        dpi=400,
        size=(640, 400),
        bottommargin=4Plots.mm,
        leftmargin=2Plots.mm,
    )
    
    
    xticks_range = [round(10. ^i, digits=abs(i)) for i in -3:2 if β_min <= 10. ^i <= β_max]
    xticks_labs = ["\$ 10^{$(i)} \$" for i in -3:2 if β_min <= 10. ^i <= β_max]
    
    yticks_range = 0:0.25:1
    yticks_labs = ["\$ $(m) \$" for m ∈ yticks_range]
    
    plot!(xticks=(xticks_range, xticks_labs),
          yticks=(yticks_range, yticks_labs),)
    
    plot!(xaxis=:log10)
    plot!(;kw...)

    return m_plot
end
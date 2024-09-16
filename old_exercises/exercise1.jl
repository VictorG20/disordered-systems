using DisorderedSystems.BeliefPropagation
using DisorderedSystems.ProjectPlots


function main()
    plot_cavity = true  # Cavity parameters vs cavity position.
    save_cavity = false  # Save the cavity parameters' plot.

    plot_magnet = true  # Mean magnetization for different β's.
    save_magnet = false  # Save mean magnetization plot.

    figs_dir = "$(splitdir(@__DIR__)[1])" * "/figures/"  # Directory to save the plots.
    
    plots = []

    """
    Notation:
    N: Number of nodes. Integer.
    β: Inverse thermal energy. Real.
    h: External magnetic field. Real.
    """

    if plot_cavity
        println("+ Plotting cavity parameters vs cavity positions...")
        N = 50 
        β = 1//2 
        h = 1//3
        probabilities::Vector{Float64} = [0.25, 0.5, 0.75, 0.9]  # First cavity prob.
        up::Bool = true   # If true, the probabilities above correspond to P(σ = +1). 
        
        # Solve for the cavity parameters and make plot.
        cavities_list = getCavityParameters.(
            N, β, h, probabilities, up=up, fill=true, print_fixed=true)
        cavity_plot = plotCavityParameters(cavities_list, probabilities)
        
        # Add to plots.
        β_label, h_label = round(β, digits=2), round(h, digits=2)
        filename = "CavityPosition_N$(N)_beta$(β_label)_h$(h_label).png"
        push!(plots, (cavity_plot, filename, save_cavity))
    end
    
    if plot_magnet
        println("+ Plotting mean magnetization at different temperatures...")
        N = 50  
        h = 1
        init_prob = 0.5  # Probability for the starting node.
        up = true  # If true, the probability above corresponds to P(σ = +1). 
        β_min, β_max = 0.01, 10.  # Minimum / maximum inverse thermal energy.
        β_length::Int = 33  # How many values (in total) between min and max.
        
        # Logarithmically equally-spaced β values.
        β_values = exp10.(range(log10(β_min), log10(β_max), length=β_length))
        m_values = getMeanMagnetization(N, β_values, h, init_prob, up)
        magnet_plot = plotMeanMagnetization(m_values, β_values, h, N)

        h_label = isinteger(h) ? h : round(h, digits=2)
        filename = "MeanMagnetization_N$(N)_h$(h_label).png"
        push!(plots, (magnet_plot, filename, save_magnet))
    end

    for (plot, filename, save_opt) in plots
        if save_opt
            figure_path = figs_dir * filename
            println("+ Saving figure to $(figure_path)")
            savefig(plot, figure_path)
        end
        display(plot)
    end

    return nothing
end

main()

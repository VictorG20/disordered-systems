using DisorderedSystems.DataHandling: viewH5FileContents
using DisorderedSystems.RandomMatrixTheory: direct_diagonalization
using DisorderedSystems.ProjectPlots: plotSpectralDensity, savefig


"""
    setPlotsKwargs()

Wrapper to set plot kwargs. These are passed as plot!(; kw...)
"""
function setPlotsKwargs()
    plots_kw = (;
        xlabel="\$ \\lambda \$",
        xminorgrid=true,
        xminorticks=2,
        xlabelfontsize=18,

        ylabel="\$ \\rho(\\lambda) \$",
        yminorgrid=true,
        yminorticks=2,
        ylabelfontsize=18,

        tickfontsize=18,
        legendtitlefontsize=16,
        legendfontsize=15,
        legend=:bottom,

        framestyle=:box,
        size=(600, 400),
        dpi=400
    )
    return plots_kw
end


"""
    getOption(plot_num::Int)

Get an array of cases based on a plot number.

A case is built as a NamedTuple with values (N, c, # of samples).
Negative values of c represent a mean connectivity of (N + c) so,
for example, c = -1 corresponds to a fully connected network.
"""
function getOption(plot_num::Int)

    cases = @NamedTuple{N::Int, c::Int, samples::Int}[]

    if plot_num == 1
        case1 = (N = 2^10, c = +3, samples = 10)
        case2 = (N = 2^10, c = -1, samples = 10)
        append!(cases, [case1, case2])
    elseif plot_num == 2
        case1 = (N = 2^12, c = -1,  samples = 1)
        case2 = (N = 2^12, c = 100, samples = 1)
        append!(cases, [case1, case2])
    else
        error("Plot number not recognized but you can create it.")
    end

    return cases
end


function main()
    # --------------------------------------------------------------------------------------
    # User options

    view_file = false  # Only view datasets. If true, nothing else is executed.

    # Option 1 corresponds to c=3 and c=N-1 for N=2^10 with 10 samples each.
    # Option 2 corresponds to c=N-1 and c=100 for N=2^12 with 1 sample each.
    # Additional options can be implemented within the function `getOption`.
    plot_num = 1

    save_figure = false
    figure_name = "direct_diago_$(plot_num).png"  # By default, saved within figures/
    plot_wigner = (plot_num == 1) ? false : true
    legendtitle = (plot_num == 1) ? "\$ N = 2^{10} \$" : "\$ N = 2^{12} \$"

    filename = "exercise2.h5" # File name to load/save data. Assumed to be within data/
    group_nm = "eigenvalues"  # Data group name within the file, can also be `nothing`.

    # --------------------------------------------------------------------------------------
    # Paths to load/save data and save figures.
    repo_dir = "$(splitdir(@__DIR__)[1])"  # Path to disordered-systems/
    data_dir = repo_dir * "/data/"  # Path to data/
    filepath = data_dir * filename

    figs_dir = repo_dir * "/figures/"  # Path to figures/
    fig_path = figs_dir * figure_name

    # --------------------------------------------------------------------------------------
    # Check the data available in the given file within the group(s) name(s).
    if view_file
        viewH5FileContents(filepath, group_names=group_nm)
        return
    end

    # --------------------------------------------------------------------------------------
    # Load or generate samples of eigenvalues for each case.
    kw = (;
        load_data = true,  # If true, look for datasets within `group_nm` in `filepath`
        save_data = true,  # If `load_data` is true and all datasets exist, this is ignored.
        get_missing = false,  # When loading data, generate missing datasets.
        filepath = filepath,  # File to either load or save data.
        group_name = group_nm,  # Group in .h5 file to either load or save data.
        seed = 42  # Random number generator seed.
    )

    cases = getOption(plot_num)
    λ_series = direct_diagonalization(cases; kw...)  # Get λ_series for the given cases.

    # --------------------------------------------------------------------------------------
    # Create histogram plot.
    kwargs = (;
        plots_kwarg = setPlotsKwargs(),  # kwargs to be passed to plot!(;kw...)
        plot_wigner = plot_wigner,  # Plot Wigner's semicircle law.
        shared_bins = true,  # If true, fix the amount of bins between min and max values.
        number_bins = 60,  # If `shared_bins`, this is the total number of bins in the plot.
        format_tick = true,  # Format axes ticks with LaTeX labels.
        legendtitle = legendtitle,
    )

    histos = []
    colors = [:blue, :forestgreen, :orange, :violet]
    histo_lw = 1.7

    for (idx, λ_samples, case, color) in zip(1:length(cases), λ_series, cases, colors)

        # A lot of words just to change a color.
        (length(cases) == 2) && !kwargs[:plot_wigner] && (idx == 2) && (color = :red)

        histo = (
            data = vec(λ_samples), 
            h_kw = (
                label = "\$ c = $(case.c < 0 ? "N - $(abs(case.c))" : case.c) \$", 
                lw=histo_lw, bins=:auto, normalize = :pdf, color=color,
                )
            )
        push!(histos, histo)
    end

    histo_plot = plotSpectralDensity(histos; kwargs...)

    if save_figure
        println("Saving figure to $(fig_path)...")
        savefig(histo_plot, fig_path)
    end
    
    display(histo_plot)

    return nothing
end


main()

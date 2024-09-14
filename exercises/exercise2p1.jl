using DisorderedSystems.DataHandling: viewH5FileContents
using DisorderedSystems.RandomMatrixTheory: direct_diagonalization
using DisorderedSystems.ProjectPlots: plotSpectralDensity


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


function main()
    view_file = false  # Only view datasets. Nothing else is executed if true.

    # Data loading/saving.
    data_dir = "$(splitdir(@__DIR__)[1])" * "/data/"  # Path to data/
    filename = "exercise2.h5"
    filepath = data_dir * filename
    group_nm = "eigenvalues"  # Data group name within the file.

    if view_file
        viewH5FileContents(filepath, group_names=group_nm)
        return
    end

    # A case is built as (N, c, # of samples).
    # Negative values of c are interpreted as mean connectivity = (N + c),
    # so, for example, c = -1 corresponds to a fully connected network.
    
    # case1 = (N = 2^10, c = +3, samples = 10)
    # case2 = (N = 2^10, c = -1, samples = 10)
    
    case1 = (N = 2^12, c = -1,  samples = 1)
    case2 = (N = 2^12, c = 100, samples = 1)
    cases::Vector{@NamedTuple{N::Int, c::Int, samples::Int}} = [case1, case2]

    kw = (;
        load_data = true,
        save_data = true,  # If `load_data` is true and all datasets exist, this is ignored.
        get_missing = false,  # When loading data, generate missing datasets.
        filepath = filepath,  # File to either load or save data.
        group_name = group_nm,  # Group in .h5 file to either load or save data.
        seed = 42  # Random number generator seed.
    )

    # Get λ_series for the given cases.
    λ_series = direct_diagonalization(cases; kw...)

    # Create histogram plot.
    histo_lw = 1.7
    histo1 = (
        data = vec(λ_series[1]), 
        h_kw = (
            color=:blue, 
            label="\$ c = $(case1.c < 0 ? "N - $(abs(case1.c))" : case1.c) \$", 
            lw=histo_lw, bins=:auto, normalize = :pdf
            )
        )
    histo2 = (
        data = vec(λ_series[2]), 
        h_kw = (
            color=:forestgreen, 
            label="\$ c = $(case2.c < 0 ? "N - $(abs(case2.c))" : case2.c) \$", 
            lw=histo_lw, bins=:auto, normalize = :pdf
            )
        )

    kwargs = (;
        plots_kwarg = setPlotsKwargs(),
        plot_wigner = true,
        shared_bins = true,
        number_bins = 60,
        format_tick = true,
        legendtitle = "\$ N = 2^{10} \$",
    )

    histo_plot = plotSpectralDensity([histo1, histo2]; kwargs...)
    
    display(histo_plot)

    return nothing
end


main()

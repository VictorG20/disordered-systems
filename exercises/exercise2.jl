using Random

using DisorderedSystems: saveData2h5, loadDataFromH5, viewH5FileContents
using DisorderedSystems.RandomMatrixTheory: sampleEigenvalues
using DisorderedSystems.ProjectPlots: plotSpectralDensity


"""
    getDefaultValues(point::Int)

Retrieve the parameter values according to each part of exercise 2 part 1.

Convenience function to avoid cluttering the main function. The points 1 to 4
correspond to the following parameter values:
1. 10 samples of linear size 2^10 with mean connectivity 3.
2. 10 samples of linear size 2^10 with mean connectivity N - 1.
3. 1 sample of linear size 2^12 with mean connectivity N - 1.
4. 1 sample of linear size 2^12 for mean connectivity values of 50, 100 and 200.
"""
function getDefaultValues(point::Int8)

    !(point in (1, 2, 3, 4)) && error("Part number not valid. Options are 1 - 4")

    if point <= 2
        num_samples = 10
        linear_size = 2^10
        avg_connect = (point == 1) ? 3 : linear_size - 1
    else
        num_samples = 1
        linear_size = 2^12
        avg_connect = (point == 3) ? linear_size - 1 : [20]
    end

    return linear_size, avg_connect, num_samples
end


function getLambdaSeries(plot_num::Int8; save_data::Bool=false, load_data::Bool=false)

    !(plot_num in [1, 2]) && error("Plot number not valid. Options are only 1 or 2.")

    points = (plot_num == 1) ? Int8[1, 2] : Int8[3, 4]

    data_dir = "$(splitdir(@__DIR__)[1])" * "/data/"  # Path to data/
    filename = "exercise2.h5"
    filepath = data_dir * filename
    group_nm = "eigenvalues"  # Data group name within the file.

    λ_series = []
    labels = []

    for point in points
        N, c_values, num_samples = getDefaultValues(point)

        for c in c_values
            data_name = "size_$(N)_c_$(c)_samples_$(num_samples)"

            if load_data
                λ_samples = loadDataFromH5(filepath, data_name; group_name=group_nm)
            else
                println("Generating $(num_samples) samples of" * 
                    " size N=$(N) and mean connectivity c=$(c)...")
                # Reseed at each value of 'c', so the order of sampling is irrelevant.
                Random.seed!(42)  
                λ_samples = sampleEigenvalues(N, c, num_samples)
                
                if save_data
                    println("  Saving data to $(filepath)")
                    saveData2h5(
                        filepath, data_name, λ_samples; 
                        group_name=group_nm, overwrite=true)
                end
            end
            push!(λ_series, λ_samples)
            push!(labels, "\$ c = $(c) \$")
        end
    end

    return λ_series, labels
end

function main()
    plot_num::Int8 = 2
    load_data = true
    save_data = true  # Ignored if values are loaded.
    # save_figure = false

    # 1.1 Generate 10 samples with N=2^10, c=3, and estimate ρ(λ).
    # 1.2 Generate 10 samples with N=2^10, c=N-1, and estimate ρ(λ). 
    # 1.3 Generate 1 sample with N=2^12, c=N-1, estimate ρ(λ) and 
    #     compare with Wigner's semicircle law.
    # 1.4 Compare spectral density histograms for different values
    #     of c, to Wigner's semicircle law.
    
    λ_series, labels = getLambdaSeries(plot_num; load_data=load_data, save_data=save_data)

    kw = (;
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

    if plot_num == 1
        push!(labels, "\$ N = $(2^10) \$")
        colors = [:blue, :red]
        plot_wigner = false
    else
        push!(labels, "\$ N = $(2^12) \$")
        colors = [:blue, :forestgreen, :red]
        plot_wigner = true
    end

    histogram = plotSpectralDensity(
        λ_series, plot_wigner=plot_wigner, colors=colors, labels=labels, kw=kw, bins=80)

    display(histogram)

    return nothing
end

# main()

direct_diagonalization()
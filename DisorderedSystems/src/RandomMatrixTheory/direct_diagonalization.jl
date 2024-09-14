"""
    direct_diagonalization(cases::Vector{@NamedTuple{N::Int, c::Int, samples::Int}}; 
    load_data::Bool=false, save_data::Bool=false, get_missing::Bool=false, 
    filepath::Union{String, Nothing}=nothing, group_name::Union{String, Nothing}=nothing, 
    seed::Union{Int, Nothing}=nothing)

Get the eigenvalues for each case in `cases`, consisting of samples of real
Gaussian Wigner Matrices of linear size ``N`` and mean connectivity ``c``


# Keyword Arguments
- `load_data::Bool=false`: Load data from `filepath` within group `group_name`.
- `save_data::Bool=false`: Save any generated data to `filepath` within group `group_name`.
- `get_missing::Bool=false`: If `load_data`, generate datasets that are not found.
- `filepath::Union{String, Nothing}=nothing`: Path to .h5 file either to save or load data.
- `group_name::Union{String, Nothing}=nothing`: Group within .h5 file to save or load data. 
    If `nothing` (default), data is loaded/saved from the top level of the file.
- `seed::Union{Int, Nothing} = nothing`: If an integer is given, the random number generator
    is seeded every time a new case starts so that the order in which data is generated is
    irrelevant.
"""
function direct_diagonalization(cases::Vector{@NamedTuple{N::Int, c::Int, samples::Int}}; 
    load_data::Bool=false, save_data::Bool=false, get_missing::Bool=false, 
    filepath::Union{String, Nothing}=nothing, group_name::Union{String, Nothing}=nothing, 
    seed::Union{Int, Nothing}=nothing)

    data_names = [
        "size_$(N)_c_$((c < 0) ? (c = N + c) : c)_samples_$(s)" for (N, c, s) in cases]

    if load_data
        println("Loading datasets...")
        λ_series = loadDataFromH5(filepath, data_names, group_name=group_name)
        miss_ids = findall(λ_series .== -1)

        # 1. All datasets existed.
        isempty(miss_ids) && (return λ_series)

        # 2. At least one dataset is missing. If not 'get_missing', throw error.
        !get_missing && error("Couldn't retrieve the following datasets:\n", 
            "  ", data_names[miss_ids], 
            "\nCheck that the data exists or set 'get_missing = true' to generate it.")
        
        # 3. At least one dataset is missing and user requested to create.
        for idx in miss_ids
            N, c, samples = cases[idx]
            (c < 0) ? (c = N + c) : c  # Correct mean connectivity from negative values.

            println("  Generating data for $(samples) samples with N=$(N) and c=$(c)...")

            # Reseed at each case so the order of the cases is irrelevant.
            (!isnothing(seed)) && Random.seed!(seed)
            λ_series[idx] = sampleEigenvalues(N, c, samples)

            if save_data
                println("    Saving data...")
                saveData2h5(filepath, data_names[idx], λ_series[idx]; 
                    group_name=group_name)
            end
        end

        return λ_series

    end  # End to load_data case.

    # If this point is reached, data is not loaded, only generated.
    λ_series = []  # Vector with each entry corresponding to the eigenvalues for each case.

    for ((N, c, samples), data_name) in zip(cases, data_names)
        (c < 0) ? (c = N + c) : c  # Correct mean connectivity from negative values.
        
        println("Generating data for $(samples) samples with N=$(N) and c=$(c)...")

        # Reseed at each case so the order of the cases is irrelevant.
        (!isnothing(seed)) && Random.seed!(seed)
        λ_samples = sampleEigenvalues(N, c, samples)
        push!(λ_series, λ_samples)

        if save_data
            println("  Saving data...")
            saveData2h5(filepath, data_name, λ_samples; group_name=group_name)
        end
    end

    return λ_series
end
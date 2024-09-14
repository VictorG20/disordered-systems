"""
    loadDataFromH5(filepath::String, data_name::String; group_name::String=nothing)

Load dataset, within a group or at the top level (default), from an .h5 file.
"""
function loadDataFromH5(filepath::String, data_name::String; 
    group_name::Union{String, Nothing}=nothing)

    fid = h5open(filepath, "r")
    group_data = isnothing(group_name) ? fid : fid[group_name] 
    dataset = read(group_data[data_name])
    close(fid)

    return dataset
end


"""
    loadDataFromH5(filepath::String, data_names::Vector{String}; 
    group_name::Union{String, Nothing}=nothing)

Load several datasets, all within the same group or all at the top level, from an .h5 file.

Returns an array of datasets with each entry corresponding to the respective data name.
If any given dataset cannot be retrieved from the file, the respective entry in the array
returns a -1.
"""
function loadDataFromH5(filepath::String, data_names::Vector{String}; 
    group_name::Union{String, Nothing}=nothing)
    
    # Open file and access group if applicable.
    fid = h5open(filepath, "r")
    group_data = isnothing(group_name) ? fid : fid[group_name]

    # Try to retrieve datasets, return -1 if key not found and throw error otherwise.
    datasets = []
    for data_name in data_names
        dataset = -1
        try
            dataset = read(group_data[data_name])
        catch error
            !isa(error, KeyError) && error("Could not retrieve $(data_name). ", error)
        end
        push!(datasets, dataset)
    end

    close(fid)

    return datasets
end
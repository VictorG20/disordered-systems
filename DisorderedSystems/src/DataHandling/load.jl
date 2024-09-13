"""
    loadDataFromH5(filepath::String, data_name::String; group_name::String=nothing)

Load data from a .h5 file.

If no group_name is given (default), the dataset is assumed to be at the top level.
"""
function loadDataFromH5(filepath::String, data_name::String; 
    group_name::Union{String, Nothing}=nothing)

    fid = h5open(filepath, "r")
    group_data = isnothing(group_name) ? fid : fid[group_name] 
    λ_samples = read(group_data[data_name])
    close(fid)

    return λ_samples
end
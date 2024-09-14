"""
    saveData2h5(filepath::String, data_name::String, data; 
        group_name::String = nothing, overwrite::Bool=false)

Save data to a dataset in an .h5 file.

The file `filepath` is created if it does not exist. Datasets are 
created at the top level of the file by default. If a group name is
given, it is created if necessary and the dataset is created within.
Setting `overwrite=true` allows to overwrite `data_name` if it already exists.
"""
function saveData2h5(filepath::String, data_name::String, data; 
    group_name::Union{String, Nothing} = nothing, overwrite::Bool=false)

    # Open file to write, creating it if it doesn't exist. 
    # Access group or create if it doesn't exist.
    fid = h5open(filepath, "cw")  

    # Access group.
    if isnothing(group_name)
        data_group = fid  # Save dataset in first level.
    elseif group_name in keys(fid)
        data_group = fid[group_name]
    else
        data_group = create_group(fid, group_name)
    end

    # Check if a dataset with the given name already exists.
    if data_name in keys(data_group)
        if ! overwrite
            error("A dataset with the name '$(data_name)' already exists. " * 
            "Try a different name or set `overwrite=true`.")
        end

        println("    Overwritting data...")
        dataset = data_group[data_name]
        delete_object(dataset)
    end

    dataset = create_dataset(data_group, data_name, Float64, size(data))
    write(dataset, data)
    close(fid)

    return nothing
end

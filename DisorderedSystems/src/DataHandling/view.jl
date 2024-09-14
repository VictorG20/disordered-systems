"""
    viewH5FileContents(filepath::String; 
    group_names::Union{String, Vector{String}, Nothing} = nothing)

View the directory tree and datasets within the .h5 file. If 
group_names are given, display only the datasets within these.
"""
function viewH5FileContents(filepath::String; 
    group_names::Union{String, Vector{String}, Nothing} = nothing)

    fid = h5open(filepath, "r")

    if isnothing(group_names)
        display(fid)
        close(fid)
        return
    end

    # If group_names is a single string, convert to vector, since Julia
    # would take the string and iterate over every character.
    isa(group_names, String) ? (group_names = [group_names]) : nothing
    
    for group_name in group_names
        try
            group = fid[group_name]
            display(group)
        catch error
            
            if isa(error, KeyError)
                println("Group name '$(group_name)' not found in file. Skipping...")
            else
                println(
                    "Something weird happened when looking for the group '$(group_name)'.", 
                    "\n    ",error
                )
            end

            continue
        end
    end

end
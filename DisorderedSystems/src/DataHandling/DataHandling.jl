module DataHandling
    using HDF5

    include("save.jl")
    export saveData2h5

    include("load.jl")
    export loadDataFromH5

    include("view.jl")
    export viewH5FileContents
    
end
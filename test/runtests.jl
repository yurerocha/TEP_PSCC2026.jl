using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Test

# Includes all test_*.jl files inside the "test" folder
function recursive_include(path::String)
    for file in readdir(path)
        file_path = joinpath(path, file)
        if isdir(file_path)
            recursive_include(file_path)
            continue
        elseif !endswith(file, ".jl")
            continue
        elseif startswith(file, "test_")
            include(file_path)
        end
    end
end

@testset verbose = true "TEP Tests" begin
    # recursive_include(@__DIR__)
    # include("test_cats.jl")
    include("test_pglibopf.jl")
    # include("test_ptdf_model.jl")
end

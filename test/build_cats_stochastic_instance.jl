using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Serialization

params = TEP.Parameters()

output_dir = "input"
file = "CaliforniaTestSystem.m"
inputfile = "external/CATS-CaliforniaTestSystem/MATPOWER/" * 
            "CaliforniaTestSystem.m"

try
    TEP.rm_dir(output_dir)
catch e
    @warn e
end

inst = TEP.build_cats_stochastic_instance(params, inputfile)
open("$output_dir/$file", "w") do io
    serialize(io, inst)
end
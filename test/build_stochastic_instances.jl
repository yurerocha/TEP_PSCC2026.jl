using TEP_PSCC2026
const TEP = TEP_PSCC2026
using Serialization

params = TEP.Parameters()

start_file = 1
end_file = 5
dir = "external/pglib-opf/"
output_dir = "input"

# try
#     TEP.rm_dir(output_dir)
# catch e
#     @warn e
# end

files = [
         "pglib_opf_case3012wp_k.m",
         "pglib_opf_case6495_rte.m",
         "pglib_opf_case7336_epigrids.m",
         "pglib_opf_case9591_goc.m",
         "pglib_opf_case10000_goc.m"
        ]

# files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
# sort!(files, by=x->parse(Int, match(r"\d+", x).match))
skip = []

for (i, file) in enumerate(files[start_file:end_file])
    inst = TEP.build_stochastic_instance(params, dir * file)
    open("$output_dir/$file", "w") do io
        serialize(io, inst)
    end
end
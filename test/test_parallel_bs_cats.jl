module TestBeamSearchCATS

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using MPI
using Random

include("utils.jl")

start_file = 1
end_file = 10
# log_dir = "test/log_bs_cats"
log_dir = "test/log"
log_file = "$log_dir/tep_bs_cats.md"
dir = "input"
base_file = "CaliforniaTestSystem"

try
    rm_dir(log_dir)
catch e
    @warn e
end

# rng = Random.MersenneTwister(123)

TEP.log_header(log_file)

files = select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))

project_dir = dirname(Base.active_project())
parallel_script = joinpath(@__DIR__, "parallel_bs_cats.jl")

for (i, file) in enumerate(files[start_file:end_file])
    mpiexec(exe -> run(`$exe -n 8 $(Base.julia_cmd()) \
                        --project=$(project_dir) $(parallel_script) \
                        $log_file $log_dir $dir $file`))
end

end # module
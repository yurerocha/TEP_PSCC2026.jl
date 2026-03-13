module TestBeamSearchPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using MPI
using Random

params = TEP.Parameters()

start_file = 40 # 40
end_file = 61 # 61
log_dir = "test/quali_ub_5min/rnf_bs"
log_file = "$log_dir/log.md"
dir = "external/pglib-opf"

try
    TEP.rm_dir(log_dir)
catch e
    @warn e
end

# rng = Random.MersenneTwister(123)

TEP.log_header(log_file)

files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))
# Run solver with binary decision variables
skip = ["pglib_opf_case8387_pegase.m"]
# "pglib_opf_case30_ieee.m" Infeasible

project_dir = dirname(Base.active_project())
parallel_script = joinpath(@__DIR__, "parallel_bs_pglibopf.jl")

for (i, file) in enumerate(files[start_file:end_file])
    if file in skip
        TEP.log(params, "Skipping instance $file")
        continue
    end
    TEP.log(params, "Test $file num $(start_file + i - 1)", true)

    mpiexec(exe -> run(`$exe -n 8 $(Base.julia_cmd()) \
                        --project=$(project_dir) $(parallel_script) \
                        $log_file $start_file $end_file $log_dir $dir $file`))
end
# mpiexec(exe -> run(`$exe -n 7 $(Base.julia_cmd()) \
#                     --project=$(project_dir) $(parallel_script)`))

end # module
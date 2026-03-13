module TestRemoveAndFixPGLibOPF

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using CSV
using DataFrames
using PowerModels
using Random
using Test

# include("utils.jl")

params = TEP.Parameters()

start_file = 40 # 40
end_file = 61 # 62
is_rnr_heur_en = true
scen = 1

log_dir = "test/quali_ub_5min/rnr"
log_file = "$log_dir/log.md"


try
    TEP.rm_dir(log_dir)
catch e
    @warn e
end

# Number of seconds since the Unix epoch
# seed = Int(floor(datetime2unix(now())))
# Random.seed!(seed)
# rng = Random.MersenneTwister(123)

dir = "external/pglib-opf"
files = TEP.select_files(dir, end_file)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))
# Run solver with binary decision variables
skip = ["pglib_opf_case8387_pegase.m"]

md_file = "test/quali_ub_5min/rnr_bs_mip1/log.md"
df = TEP.read_markdown_table(md_file)
# Solve (s): 8
# RNR (s): 15
# BS (s): 19
time_limit = parse.(Float64, df[!, 8]) + 
             parse.(Float64, df[!, 16]) + 
             parse.(Float64, df[!, 19])
# time_limit = []

TEP.log_header(log_file)

@testset "[RNR Heur] PGLibOPF" begin
    for (i, file) in enumerate(files[start_file:end_file])
        if file in skip
            TEP.log(params, "Skipping instance $file")
            continue
        end
        tl = max(time_limit[21], params.solver.time_limit)
        # tl = params.solver.time_limit
        TEP.log(params, "Processing $file $(start_file + i - 1) $tl", true)

        filepath = "$(dir)/$file"
        logsolver = "$(dir)/log/$file"

        params.log_dir = log_dir
        params.log_file = "$log_dir/$file"
        params.solver.time_limit = tl

        inst = TEP.build_instance(params, filepath)

        mip = nothing
        build_time = 0.0
        start_time = 0.0
        results = TEP.init_results()

        try
            TEP.log(params, "Build full model", true)
            build_time = @elapsed (mip = TEP.build_mip(inst, params))

            lp = TEP.build_lp(inst, params, scen)
            TEP.log(params, "Build heuristic solution", true)
            (start, report) = TEP.build_solution!(inst, params, scen, 
                                                  lp, false, TEP.Cache(0, 0), 
                                                  is_rnr_heur_en)

            TEP.log(params, "Fix the start of the model", true)
            fix_start_time = 
                @elapsed (start_ub = TEP.fix_start!(inst, params, mip, start))

            params.solver.time_limit -= (fix_start_time + report.time)

            # # TODO: Add tuning flag
            # TEP.log(params, "Solve the model", true)
            # results = TEP.solve!(inst, params, mip)

            # results["build_time"] = build_time
            # if isa(results["objective"], Number)
            #     results["build_obj_rat"] = TEP.comp_build_obj_rat(inst, 
            #                                                 results["objective"], 
            #                                                 start.inserted)
            # end

            results["rnr_time"] = report.time
            results["rnr_rm_rat"] = report.rm_ratio
            results["rnr_impr_rat"] = report.impr_ratio
            results["fix_start_time"] = fix_start_time
            results["start_ub"] = start_ub
            
            TEP.log_instance(log_file, file, inst, results)
        catch e
            @warn e
            TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
        end
    end
end # testset

end # module
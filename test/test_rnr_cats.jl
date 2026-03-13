module TestRemoveAndFixCATS

using TEP_PSCC2026
const TEP = TEP_PSCC2026
using CSV
using DataFrames
using PowerModels
using Random
using Test

include("utils.jl")

params = TEP.Parameters()
params.instance.cost_mult = 1e3

start_file = 1
end_file = 10
num_scenarios = 10
dir = "input"
base_file = "CaliforniaTestSystem"
is_rnr_heur_en = false

suffix = ""
if is_rnr_heur_en
    suffix = "_rnr"
end
log_dir = "test/log$(suffix)_cats"
log_file = "$log_dir/tep$(suffix)_cats.md"

try
    rm_dir(log_dir)
catch e
    @warn e
end

# rng = Random.MersenneTwister(123)

dir = "input"
files = select_files(dir, num_scenarios)
# Sort files so that the smallest instances are solved first
sort!(files, by=x->parse(Int, match(r"\d+", x).match))

md_file = "test/log_bs_cats/tep_bs_cats.md"
df = read_markdown_table(md_file)
# Solve (s): 8
# RNR (s): 15
# BS (s): 19
time_limit = parse.(Float64, df[!, 8]) + 
             parse.(Float64, df[!, 15]) + 
             parse.(Float64, df[!, 19])

TEP.log_header(log_file)

@testset "[RNR Heur] CATS" begin
    for (i, file) in enumerate(files[start_file:end_file])
        tl = time_limit[i]
        TEP.log(params, "Processing $file $(start_file + i - 1) $tl", true)
    
        filepath = "$(dir)/$file"
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
    
            TEP.log(params, "Build heuristic solution", true)
            (start, report) = TEP.build_solution(inst, params, is_rnr_heur_en)
    
            TEP.log(params, "Fix the start of the model", true)
            fix_start_time = 
                    @elapsed (obj = TEP.fix_start!(inst, params, mip, start))
            results["objective"] = obj
    
            params.solver.time_limit -= (fix_start_time + report.time)
    
            # TODO: Add tuning flag
            TEP.log(params, "Solve the model", true)
            results = TEP.solve!(inst, params, mip)
    
            results["build_time"] = build_time
            if isa(results["objective"], Number)
                results["build_obj_rat"] = TEP.comp_build_obj_rat(inst, 
                                                            results["objective"], 
                                                            start.inserted)
            end
    
            results["rnr_time"] = report.time
            results["rnr_rm_percent"] = report.removed_percent
            results["rnr_impr_percent"] = report.improvement_percent
            results["fix_start_time"] = fix_start_time
            
            TEP.log_instance(log_file, file, inst, results)
        catch e
            @warn e
            TEP.log_instance(log_file, "<s>" * file * "</s>", inst, Dict())
        end
    end
end # testset

end # module
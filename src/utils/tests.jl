"""
    select_files(path::String, num_files::Int64)

Benchmark: https://github.com/power-grid-lib/pglib-opf
"""
function select_files(path::String, num_files::Int64)
    files = []
    for file in readdir(path)
        if endswith(file, ".m")
            push!(files, file)
        end
    end
    # Sort files so that instances with less buses are sovled first
    sort!(files, by=x->parse(Int, match(r"\d+", x).match))
    return files[1:num_files]
end

function read_markdown_table(file_path::String)
    lines = readlines(file_path)

    # Filter out lines that are not part of the table (remove the separator row with dashes)
    table_lines = filter(line -> occursin('|', line) && 
                                !occursin(":---", line), lines)

    # Strip outer pipes and split into fields
    rows = [split(line, "|")[2:end - 1] for line in table_lines]

    # Create a DataFrame
    header = strip.(rows[1])
    data = [strip.(r) for r in rows[2:end]]

    return DataFrame([Symbol(h) => 
                      [row[i] for row in data] for (i, h) in enumerate(header)])
end

function rm_dir(dir::String)
    if isdir(dir)
        if !isempty(dir)
            rm(dir, recursive = true)
        end
    end
    mkdir(dir)
end
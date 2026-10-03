"""
    PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               write_log = false,
                               save_diags = false,
                               max_filesize = Inf,
                               )
Generate a `PlanktonOutputWriter` structure which includes settings for model outputs

Keyword Arguments (Optional)
============================
- `dir`: The directory to store model outputs, "./results" by default
- `diags_prefix`: Descriptive filename prefixed to diagnostic output files.
- `write_log`: write model logs which contain global averages of simulated phytoplankton, default: `false`.
- `save_diags`: write diagnostics to disk, default: `false`.
- `max_filesize`: The writer will stop writing to the output file once the file size exceeds `max_filesize`,
                    and write to a new one with a consistent naming scheme ending in `part1`, `part2`, etc.
                    default: `Inf`.
"""
function PlanktonOutputWriter(;dir = "./results",
                               diags_prefix = "diags",
                               write_log = false,
                               save_diags = false,
                               max_filesize = Inf
                               )

    max_filesize > 0 || throw(ArgumentError("max_filesize must be positive"))
    isdir(dir) && rm(dir, recursive=true)
    mkpath(dir)

    diags_file = ""

    if save_diags
        diags_file = joinpath(dir, diags_prefix*".jld2")
    end

    return PlanktonOutputWriter(dir, write_log, save_diags, diags_file, max_filesize, 1)
end

function show(io::IO, writer::PlanktonOutputWriter)
    print(io, "PlanktonOutputWriter:\n",
              "├── files are saved at $(writer.filepath)\n",
              "├── $(save_diags_string(writer))\n",
              "├── write log: $(writer.write_log)\n",
              "└── Maximum file size: $(humanize_filesize(writer.max_filesize))"
              )
end

function save_diags_string(writer::PlanktonOutputWriter)
    if writer.save_diags
        return "diagnostics are saved as $(writer.diags_file)"
    else
        return "diagnostics are not saved"
    end
end

function humanize_filesize(s::Number)
    isinf(s) && return "unlimited"
    suffix = ["B", "KiB", "MiB", "GiB", "TiB", "PiB", "EiB", "ZiB", "YiB"]
    biggest_suffix = suffix[1]
    value = s
    unit = 1024.0f0
    for power in eachindex(suffix)
        unit = 1024.0^power
        biggest_suffix = suffix[power]
        value < unit && break
    end
    s = 1024 * value / unit
    return @sprintf("%3.1f %s", s, biggest_suffix)
end

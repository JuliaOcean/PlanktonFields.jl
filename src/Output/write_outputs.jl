##### write model outputs
function write_output!(writer::Union{PlanktonOutputWriter, Nothing}, model::PlanktonModel, diags::PlanktonDiagnostics)
    if isa(writer, Nothing)
        return nothing
    else
        if writer.write_log
            write_species_dynamics(model.t, model.plankton.phytos,
                                   writer.filepath, model.mode)
        end

        if writer.save_diags
            if model.iteration % diags.iteration_interval == 0.0f0
                if isfile(writer.diags_file) && filesize(writer.diags_file) ≥ writer.max_filesize
                    start_next_diags_file(writer)
                end
                write_diags_to_jld2(diags, writer.diags_file, model.t, model.iteration,
                                    diags.iteration_interval, model.grid)
            end
        end
    end
end

function start_next_diags_file(writer::PlanktonOutputWriter)
    if writer.part_diags == 1
        part1_path = replace(writer.diags_file, r".jld2$" => "_part1.jld2")
        mv(writer.diags_file, part1_path)
        writer.diags_file = part1_path
    end

    writer.part_diags += 1
    writer.diags_file = replace(writer.diags_file, r"part\d+.jld2$" => "part" * string(writer.part_diags) * ".jld2")
end

function write_diags_to_jld2(diags, filepath, t, iter, ncounts, grid)
    jldopen(filepath, "a+") do file
        file["timeseries/t/$iter"] = t
        for key in keys(diags.bgc_tracers)
            file["timeseries/$key/$iter"] = Array(interior(diags.bgc_tracers[key], grid)) ./ ncounts
        end
        for (group, species) in (("phytoplankton", diags.phyto_tracers), ("zooplankton", diags.zoo_tracers))
            for sp in keys(species), proc in keys(species[sp])
                file["timeseries/$group/$sp/$proc/$iter"] = Array(interior(species[sp][proc], grid)) ./ ncounts
            end
        end
    end
    ##### zeros diags
    for tr in diags.bgc_tracers
        tr .= 0.0f0
    end
    for species in (diags.phyto_tracers, diags.zoo_tracers), sp in species, proc in sp
        proc .= 0
    end
    return nothing
end

function write_species_dynamics(args...)
    throw(ArgumentError("Species logging is not implemented yet"))
end

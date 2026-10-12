"""
    TimeStep!(model::PlanktonModel, ΔT::Int64, diags::PlanktonDiagnostics, resultpath::String)
Update physiology processes and tracer field of `PlanktonModel` one time step forward.

Keyword Arguments
=================
- `model`: `PlanktonModel` to be updated one time step forward.
- `ΔT`: The length of a time step.
- `diags`: `PlanktonDiagnostics` to be updated.
"""
function TimeStep!(model::PlanktonModel, ΔT, diags::PlanktonDiagnostics)
    model.t = model.t+ΔT
    model.iteration = model.iteration+1
 
    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.Chl .= 0.0f0

    ##### calculate total chlorophyll for PAR attenuation
    for sp in eachindex(model.plankton.phytos)
        model.timestepper.Chl .+= model.plankton.phytos[sp].fields.Chl.data
    end
    ##### calculate PAR
    for ki in 1:model.grid.Nz
        calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, 
                  model.timestepper.PARF, model.grid, model.options.kc,
                  model.options.kw, ki)
    end # PAR

    ##### update plankton fields
    plankton_update!(model.plankton, model.bgc_tracers, 
                     model.timestepper, model.grid, ΔT,
                     model.arch, model.mode)

    ##### bgc tracers update
    bgc_tracer_update!(model.bgc_tracers, model.timestepper.bgc_Gcs, 
                       model.timestepper.bgc_tmp,
                       model.timestepper.f_sink, model.arch,
                       model.grid, model.bgc_params, model.timestepper.vel, model.timestepper.plk, ΔT,
                       model.iteration)

    ##### diagnostics for plankton fields
    for sp in eachindex(model.plankton.phytos)
        for key in eachindex(diags.phytos[sp])
            if key in eachindex(model.plankton.phytos[sp].fields)
                @inbounds diags.phytos[sp][key] .+= model.plankton.phytos[sp].fields[key].data
            end
            if key in eachindex(model.plankton.phytos[sp].procs)
                @inbounds diags.phytos[sp][key] .+= model.plankton.phytos[sp].procs[key].data
            end
        end
    end
    for sp in eachindex(model.plankton.zoos)
        for key in eachindex(diags.zoos[sp])
            if key in eachindex(model.plankton.zoos[sp].fields)
                @inbounds diags.zoos[sp][key] .+= model.plankton.zoos[sp].fields[key].data
            end
            if key in eachindex(model.plankton.zoos[sp].procs)
                @inbounds diags.zoos[sp][key] .+= model.plankton.zoos[sp].procs[key].data
            end
        end
    end
    
    ##### diagnostics for bgc tracers
    @inbounds diags.bgc_tracers.PAR .+= model.timestepper.par
    for key in eachindex(diags.bgc_tracers)
        if key in eachindex(model.bgc_tracers)
            @inbounds diags.bgc_tracers[key] .+= model.bgc_tracers[key].data
        end
    end # tracers

    @inbounds model.timestepper.par₀ .= model.timestepper.par

    return nothing
end
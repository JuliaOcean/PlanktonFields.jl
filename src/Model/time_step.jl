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
 

    @inbounds model.timestepper.vel½.u.data .= (model.timestepper.vel₀.u.data .+ model.timestepper.vel₁.u.data) .* 0.5f0
    @inbounds model.timestepper.vel½.v.data .= (model.timestepper.vel₀.v.data .+ model.timestepper.vel₁.v.data) .* 0.5f0
    @inbounds model.timestepper.vel½.w.data .= (model.timestepper.vel₀.w.data .+ model.timestepper.vel₁.w.data) .* 0.5f0

    zero_fields!(model.timestepper.plk)
    @inbounds model.timestepper.Chl .= 0.0f0

    # calculate total chlorophyll for PAR attenuation
    for sp in eachindex(model.plankton.phytos)
        model.timestepper.Chl .+= model.plankton.phytos[sp].data.Chl.data
    end

    ##### calculate PAR
    for ki in 1:model.grid.Nz
        calc_par!(model.timestepper.par, model.arch, model.timestepper.Chl, 
                  model.timestepper.PARF, model.grid, model.options.kc,
                  model.options.kw, ki)
    end # PAR

    for sp in eachindex(model.plankton.phytos)
        phyto_update!(model.plankton.phytos[sp], model.bgc_tracers, 
                      model.timestepper.plk, model.timestepper.par, 
                      model.timestepper.temp, model.arch, model.grid, 
                      ΔT, model.mode)
    end

    ##### diagnostics for plankton tracers
    for sp in eachindex(model.plankton.phytos)
        for proc in keys(model.plankton.phytos[sp].data)
            @inbounds diags.phyto_tracers[sp][proc] .+= model.plankton.phytos[sp].data[proc].data
        end
    end
    for sp in eachindex(model.plankton.zoos)
        for proc in keys(model.plankton.zoos[sp].data)
            @inbounds diags.zoo_tracers[sp][proc] .+= model.plankton.zoos[sp].data[proc].data
        end
    end
    
    ##### tracers update
    bgc_tracer_update!(model.bgc_tracers, model.timestepper.bgc_Gcs, model.timestepper.bgc_tracer_tmp,
                       model.timestepper.flux_sink, model.arch,
                       model.grid, model.bgc_params, model.timestepper.vel₁, model.timestepper.plk, ΔT,
                       model.iteration)

    ##### diagnostics for bgc tracers
    @inbounds diags.bgc_tracers.PAR .+= model.timestepper.par
    for key in eachindex(diags.bgc_tracers)
        if key in eachindex(model.bgc_tracers)
            @inbounds diags.bgc_tracers[key] .+= model.bgc_tracers[key].data
        end
    end # tracers

    @inbounds model.timestepper.vel₀.u.data .= model.timestepper.vel₁.u.data
    @inbounds model.timestepper.vel₀.v.data .= model.timestepper.vel₁.v.data
    @inbounds model.timestepper.vel₀.w.data .= model.timestepper.vel₁.w.data
    @inbounds model.timestepper.par₀ .= model.timestepper.par

    return nothing
end

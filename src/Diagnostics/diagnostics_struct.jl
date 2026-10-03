"""
    PlanktonDiagnostics(model; tracer=(:PAR, :NH4, :NO3, :DOC),
                               phytoplankton = (),
                               zooplankton = (),
                               time_interval = 1)

Generate a `PlanktonDiagnostics` structure.

Keyword Arguments (Optional)
============================
- `tracer` : a `Tuple` containing the names of nutrient fields to be diagnosed.
- `phytoplankton` : a `Tuple` containing the names of physiological processes of phytoplankton individuals to be diagnosed.
- `zooplankton` : a `Tuple` containing the names of physiological processes of zooplankton individuals to be diagnosed.
- `iteration_interval` : The number of timesteps that diagnostics is averaged, 1 iteration by default.
"""
function PlanktonDiagnostics(model; tracer=(),
                                    phytoplankton = (),
                                    zooplankton = (),
                                    iteration_interval::Int = 1)
    
    iteration_interval > 0 || throw(ArgumentError("Diagnostic interval must be positive"))
    length(unique(tracer)) == length(tracer) || throw(ArgumentError("Duplicate tracer diagnostics"))
    @assert isa(tracer, Tuple)
    @assert isa(phytoplankton, Tuple)
    @assert isa(zooplankton, Tuple)

    diag_avail(tracer, phytoplankton, zooplankton, model)

    ntr   = length(tracer)
    nproc_phyto = length(phytoplankton)
    nproc_zoo = length(zooplankton)

    trs   = []
    phyto_procs = []
    zoo_procs = []
    FT = model.FT
    total_size = (model.grid.Nx+model.grid.Hx*2, model.grid.Ny+model.grid.Hy*2, model.grid.Nz+model.grid.Hz*2)

    ##### tracers
    for i in 1:ntr
        tr = zeros(FT, total_size) |> array_type(model.arch)
        push!(trs, tr)
    end
    tr_d1 = zeros(FT, total_size) |> array_type(model.arch)
    tr_default = (PAR = tr_d1,)
    diag_tr = NamedTuple{tracer}(trs)
    diag_tr = merge(diag_tr, tr_default) # add PAR as default diagnostic

    ##### phytoplankton
    phyto_name = keys(model.plankton.phytos)
    Nsp_phyto = length(phyto_name)
    for i in 1:Nsp_phyto
        procs_sp = []
        for j in 1:nproc_phyto
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{phytoplankton}(procs_sp)
        push!(phyto_procs, diag_proc)
    end
    diag_phyto = NamedTuple{phyto_name}(phyto_procs)

    ##### zooplankton
    zoo_name = keys(model.plankton.zoos)
    Nsp_zoo = length(zoo_name)
    for i in 1:Nsp_zoo
        procs_sp = []
        for j in 1:nproc_zoo
            proc = zeros(FT, total_size) |> array_type(model.arch)
            push!(procs_sp, proc)
        end
        diag_proc = NamedTuple{zooplankton}(procs_sp)
        push!(zoo_procs, diag_proc)
    end
    diag_zoo = NamedTuple{zoo_name}(zoo_procs)

    diagnostics = PlanktonDiagnostics(diag_phyto, diag_zoo, diag_tr, iteration_interval)

    return diagnostics
end

function show(io::IO, diags::PlanktonDiagnostics)
    print(io, "PlanktonDiagnostics:\n",
              "├── diagnostics of BGC tracers: $(keys(diags.bgc_tracers))\n",
              "├── diagnostics of phytoplankton: $(keys(diags.phyto_tracers))\n",
              "├── diagnostics of zooplankton: $(keys(diags.zoo_tracers))\n",
              "└── save averaged diagnostics every $(diags.iteration_interval) timesteps")
end

function diag_avail(tracer, plank, zoo, model)
    tracer_avail = tracer_avail_diags()
    plank_avail = keys(model.plankton.phytos.sp1.data)
    zoo_avail = keys(model.plankton.zoos.sp1.data)

    for i in eachindex(tracer)
        if tracer[i] ∉ tracer_avail
            throw(ArgumentError("$(tracer[i]) is not one of the diagnostics of bgc tracers"))
        end
    end

    for i in eachindex(plank)
        if plank[i] ∉ plank_avail
            throw(ArgumentError("$(plank[i]) is not one of the diagnostics of phytoplankton"))
        end
    end

    for i in eachindex(zoo)
        if zoo[i] ∉ zoo_avail
            throw(ArgumentError("$(zoo[i]) is not one of the diagnostics of zooplankton"))
        end
    end
end

function tracer_avail_diags()
    return (:PAR, :DIC, :DOC, :POC, :NH4, :NO3, :O2, :DON, :PON, :PO4, :DOP, :POP, :DFe, :PIFe, :POFe, :Dust)
end

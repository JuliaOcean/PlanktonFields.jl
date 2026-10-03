function construct_phytoplankton(arch::Architecture, sp::Int, params::Dict, FT::DataType, g::AbstractGrid)
    # construct phytoplankton tracer fields
    phyto_tracer_names = (:Chl, :C, :PS, :RS, :mort, :graz)
    data = init_tracers(arch, g, phyto_tracer_names, FT)

    param_names=(:Cquota, :Chl2C, :α, :Φ, :Topt, :Tmax, :Ea_ref, :Topt_ref, :PCmax, :respir, :k_mort, :grazFracC, :mortFracC)

    # construct phytoplankton parameters
    pkeys = Symbol.(collect(keys(params)))
    tmp = zeros(length(param_names))
    for i in eachindex(param_names)
        if param_names[i] ∉ pkeys
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(FT.(tmp))

    return PlanktonTracer(data, p)
end

function construct_zooplankton(arch::Architecture, sp::Int, params::Dict, FT::DataType, g::AbstractGrid)
    # construct zooplankton tracer fields
    zoo_tracer_names = (:C, :G, :RS, :mort,)
    data = init_tracers(arch, g, zoo_tracer_names, FT)

    param_names=(:Cquota, :Topt, :Tmax, :Ea_ref, :Topt_ref, :Gmax, :respir, :k_mort,  :mortFracC)

    # construct zooplankton parameters
    pkeys = Symbol.(collect(keys(params)))
    tmp = zeros(length(param_names))
    for i in eachindex(param_names)
        if param_names[i] ∉ pkeys
            throw(ArgumentError("PARAM: parameter not found $(param_names[i])"))
        else
            tmp[i] = params[string(param_names[i])][sp]
        end
    end
    p = NamedTuple{param_names}(FT.(tmp))

    return PlanktonTracer(data, p)
end

function initialize_phytoplankton!(phytos, source::NamedTuple, grid::AbstractGrid, arch::Architecture, FT::DataType)
    phyto_names = keys(phytos)
    for sp in eachindex(phyto_names)
        if !(phyto_names[sp] ∈ keys(source))
            throw(ArgumentError("Source for phytoplankton $(phyto_names[sp]) not found"))
        end
        phyto = phytos[sp].data
        fieldnames =  (:Chl, :C)
        for f in fieldnames
            if !(f ∈ keys(source[sp]))
                throw(ArgumentError("Source for phytoplankton $(phyto_names[sp]) must contain $(f)"))
            end
            initialize_tracer!(phyto[f], grid, source[sp][f], arch, FT)
            @views @. phyto[f].data *= grid.landmask
        end
        fill_halo_tracer!(phyto, grid)
    end
end

function initialize_zooplankton!(zoos, source::NamedTuple, grid::AbstractGrid, arch::Architecture, FT::DataType)
    zoo_names = keys(zoos)
    for sp in eachindex(zoo_names)
        if !(zoo_names[sp] ∈ keys(source))
            throw(ArgumentError("Source for zooplankton $(zoo_names[sp]) not found"))
        end
        zoo = zoos[sp].data
        fieldnames =  (:C,)
        for f in fieldnames
            if !(f ∈ keys(source[sp]))
                throw(ArgumentError("Source for zooplankton $(zoo_names[sp]) must contain $(f)"))
            end
            initialize_tracer!(zoo[f], grid, source[sp][f], arch, FT)
            @views @. zoo[f].data *= grid.landmask
        end
        fill_halo_tracer!(zoo, grid)
    end
end

function phyto_tracer_default_init()
    Chl1 = (init = 0.1, rand_noise = 0.1)
    C1 = (init = 1.0, rand_noise = 0.1)
    sp1 = (Chl = Chl1, C = C1)
    return (sp1 = sp1,)
end

function zoo_tracer_default_init()
    C1 = (init = 1.0, rand_noise = 0.1)
    sp1 = (C = C1,)
    return (sp1 = sp1,)
end
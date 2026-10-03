module Plankton

export phyto_params_default, update_phyto_params, zoo_params_default, update_zoo_params

using KernelAbstractions

using PlanktonKernels.Architectures: device, Architecture
using PlanktonKernels.Fields: Field, fill_halo_tracer!
using PlanktonKernels.Grids: AbstractGrid, volume, ΔzF

using PlanktonFields: CarbonMode, AbstractMode, PlanktonTracer, PlanktonTracers, plankton_setup

include("plankton_params.jl")
include("CarbonMode/CarbonMode.jl")
include("utils.jl")

import .Carbon

#####
##### generate plankton of multiple species
#####
function generate_plankton_tracers(plankton::plankton_setup, arch::Architecture, FT::DataType, grid::AbstractGrid, mode::AbstractMode)
    phyto_names = []
    phytos =[]
    for sp in 1:plankton.phyto_Nsp
        name = Symbol("sp"*string(sp))
        phyto = construct_phytoplankton(arch, sp, plankton.phyto_params, FT, grid, mode)
        push!(phyto_names, name)
        push!(phytos, phyto)
    end
    phytos = NamedTuple{Tuple(phyto_names)}(phytos)
    initialize_phytoplankton!(phytos, plankton.phyto_init, grid, arch, FT, mode)

    zoo_names = []
    zoos = []
    for sp in 1:plankton.zoo_Nsp
        name = Symbol("sp"*string(sp))
        zoo = construct_zooplankton(arch, sp, plankton.zoo_params, FT, grid, mode)
        push!(zoo_names, name)
        push!(zoos, zoo)
    end
    zoos = NamedTuple{Tuple(zoo_names)}(zoos)
    initialize_zooplankton!(zoos, plankton.zoo_init, grid, arch, FT, mode)

    return PlanktonTracers(phytos, zoos)
end

#####
##### some workarounds for function names
#####
construct_phytoplankton(arch::Architecture, sp::Int, params::Dict, FT::DataType, g::AbstractGrid, mode::CarbonMode) = Carbon.construct_phytoplankton(arch, sp, params, FT, g)
construct_zooplankton(arch::Architecture, sp::Int, params::Dict, FT::DataType, g::AbstractGrid, mode::CarbonMode) = Carbon.construct_zooplankton(arch, sp, params, FT, g)

initialize_phytoplankton!(phytos, source::NamedTuple, grid::AbstractGrid, arch::Architecture, FT::DataType, mode::CarbonMode) = Carbon.initialize_phytoplankton!(phytos, source, grid, arch, FT)
initialize_zooplankton!(zoos, source::NamedTuple, grid::AbstractGrid, arch::Architecture, FT::DataType, mode::CarbonMode) = Carbon.initialize_zooplankton!(zoos, source, grid, arch, FT)

phyto_tracer_default_init(mode::CarbonMode) = Carbon.phyto_tracer_default_init()
zoo_tracer_default_init(mode::CarbonMode) = Carbon.zoo_tracer_default_init()

phyto_update!(phyto::PlanktonTracer, bgc_tracers::NamedTuple, plk, par, temp, arch::Architecture, grid::AbstractGrid, ΔT::Real, mode::CarbonMode) = Carbon.phyto_update!(phyto, bgc_tracers, plk, par, temp, arch, grid, ΔT, mode)

end
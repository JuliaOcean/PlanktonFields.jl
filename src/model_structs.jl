"""
    AbstractMode
Abstract type for phytoplankton physiology modes supported by PlanktonFields.
"""
abstract type AbstractMode end

"""
    CarbonMode <: AbstractMode
Type for the phytoplankton physiology mode which only resolves carbon quota.
"""
struct CarbonMode <: AbstractMode end

"""
    QuotaMode <: AbstractMode
Type for the phytoplankton physiology mode which resolves carbon, nitrogen, and phosphorus quotas.
"""
struct QuotaMode <: AbstractMode end

##### struct that holds all plankton species
mutable struct PlanktonTracers
    phytos::NamedTuple
    zoos::NamedTuple
end
##### struct for plankton (phytoplankton/zooplankton)
##### similar to Field in PlanktonKernels, but bc replaced with a NamedTuple of parameters
mutable struct PlanktonTracer
    data::NamedTuple
    p::NamedTuple
end

mutable struct plankton_setup
    phyto_params::Union{Nothing, Dict}
    zoo_params::Union{Nothing, Dict}
    phyto_init::NamedTuple
    zoo_init::NamedTuple
    phyto_Nsp::Int64
    zoo_Nsp::Int64
end

mutable struct timestepper
    bgc_Gcs::NamedTuple                         # a NamedTuple same as bgc tracers to store tendencies
    bgc_tracer_tmp::NamedTuple                  # a NamedTuple same as bgc tracers to store tracers fields in multi-dims advection scheme
    plankton_Gcs::NamedTuple                    # a NamedTuple same as plankton tracers to store tendencies
    plankton_tracer_tmp::NamedTuple             # a NamedTuple same as plankton tracers to store tracers fields in multi-dims advection scheme
    vel₀::NamedTuple                            # a NamedTuple with u, v, w velocities
    vel½::NamedTuple                            # a NamedTuple with u, v, w velocities
    vel₁::NamedTuple                            # a NamedTuple with u, v, w velocities
    PARF::AbstractArray                         # a (Cu)Array to store surface PAR field of each timestep
    temp::AbstractArray                         # a (Cu)Array to store temperature field of each timestep
    flux_sink::AbstractArray                    # a (Cu)Array to store sinking flux field of each timestep
    plk::NamedTuple                             # a NamedTuple same as bgc tracers to store interactions with plankton fields
    par::AbstractArray                          # a (Cu)Array to store PAR field of each timestep
    par₀::AbstractArray                         # a (Cu)Array to store PAR field of the previous timestep
    Chl::AbstractArray                          # a (Cu)Array to store total chlorophyll field of each timestep
end

mutable struct ModelOpts
    kc::Float64                 # light attenuation coefficient, unit: m^-1
    kw::Float64                 # light attenuation coefficient, unit: m^-1
end

mutable struct PlanktonModel
    arch::Architecture                   # architecture on which models will run
    options::ModelOpts                   # model options
    FT::DataType                         # floating point data type
    t::AbstractFloat                     # time in second
    iteration::Int                       # model interation
    plankton::PlanktonTracers            # Tracer fields for plankton
    bgc_tracers::NamedTuple              # BGC tracer fields
    grid::AbstractGrid                   # grid information
    bgc_params::Dict                     # biogeochemical parameter set
    timestepper::timestepper             # operating Tuples and arrays for timestep
    mode::AbstractMode                   # Carbon, Quota, or MacroMolecular
end

mutable struct PlanktonDiagnostics
    phyto_tracers::NamedTuple          # for each species of phytoplankton
    zoo_tracers::NamedTuple            # for each species of zooplankton
    bgc_tracers::NamedTuple            # for tracers
    iteration_interval::Int            # time interval that the diagnostics is time averaged
end

mutable struct PlanktonOutputWriter
    filepath::String
    write_log::Bool
    save_diags::Bool
    diags_file::String
    max_filesize::Number # in Bytes
    part_diags::Int
end

mutable struct PlanktonInput
    temp::AbstractArray{<:AbstractFloat,4}      # temperature
    PARF::AbstractArray{<:AbstractFloat,3}      # PARF
    vels::NamedTuple                          # velocity fields
    ΔT_vel::AbstractFloat                     # time step of velocities provided
    ΔT_PAR::AbstractFloat                     # time step of surface PAR provided
    ΔT_temp::AbstractFloat                    # time step of temperature provided
end

mutable struct PlanktonSimulation
    model::PlanktonModel                                # Model object
    input::PlanktonInput                                # model input, temp, PAR, and velocities
    diags::Union{PlanktonDiagnostics,Nothing}           # diagnostics
    ΔT::AbstractFloat                                   # model time step
    iterations::Int                                     # run the simulation for this number of iterations
    output_writer::Union{PlanktonOutputWriter,Nothing}  # Output writer
end

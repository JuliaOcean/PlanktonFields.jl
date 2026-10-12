module Model

export model_options

import PlanktonFields: PlanktonModel, timestepper
using PlanktonFields: AbstractMode, CarbonMode, plankton_setup, ModelOpts, PlanktonDiagnostics

using PlanktonKernels.Architectures: CPU, Architecture, isfunctional, array_type
using PlanktonKernels.Grids: AbstractGrid, short_show, replace_grid_storage
using PlanktonKernels.Fields: Field, init_tracers, zero_fields!
using PlanktonKernels.Biogeochemistry: bgc_tracer_names, bgc_tracer_default_init, generate_bgc_tracers, update_bgc_params, bgc_tracer_update!

using PlanktonFields.Plankton: generate_plankton_tracers, phyto_tracer_default_init, zoo_tracer_default_init, phyto_params_default, update_phyto_params, zoo_params_default, update_zoo_params, plankton_update!, calc_par!

import Base: show

include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end

module Model

import PlanktonFields: AbstractMode, CarbonMode, QuotaMode, plankton_setup, PlanktonModel, ModelOpts, timestepper, PlanktonDiagnostics

export PlanktonModel, model_options

using PlanktonKernels.Architectures: CPU, GPU, Architecture, isfunctional, array_type
using PlanktonKernels.Grids: AbstractGrid, short_show, replace_grid_storage
using PlanktonKernels.Fields: Field, init_tracers, fill_halo_tracer!, zero_fields!
using PlanktonKernels.Biogeochemistry: bgc_tracer_names, bgc_tracer_default_init, generate_bgc_tracers, update_bgc_params, bgc_tracer_update!

using PlanktonFields.Plankton: generate_plankton_tracers, generate_plankton_tracers, phyto_tracer_default_init, zoo_tracer_default_init, phyto_params_default, update_phyto_params, zoo_params_default, update_zoo_params, phyto_update!, calc_par!

import Base: show

include("timestepper.jl")
include("models.jl")
include("time_step.jl")

end

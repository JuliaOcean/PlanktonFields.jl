module Simulation

export update!, set_vels_fields!, set_PARF_fields!, set_temp_fields!, default_PARF

import PlanktonFields: PlanktonInput, PlanktonSimulation

using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids: AbstractGrid, Bounded, replace_grid_storage
using PlanktonKernels.Fields: validate_bcs, vel_copy!, copy_interior!

using PlanktonFields: PlanktonModel, PlanktonDiagnostics, PlanktonOutputWriter
using PlanktonFields.Output: humanize_filesize, write_output!
using PlanktonFields.Model: TimeStep!

import Base: show

include("simulations.jl")
include("utils.jl")
include("update.jl")

end

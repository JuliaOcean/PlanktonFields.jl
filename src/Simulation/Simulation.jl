module Simulation
import PlanktonFields: PlanktonInput, PlanktonSimulation
using PlanktonFields: PlanktonModel, PlanktonDiagnostics, PlanktonOutputWriter
using PlanktonKernels.Architectures: CPU
using PlanktonKernels.Grids: AbstractGrid, Bounded, replace_grid_storage
using PlanktonKernels.Fields: validate_bcs, vel_copy!, copy_interior!
using PlanktonFields.Output: humanize_filesize, write_output!
using PlanktonFields.Model: TimeStep!
import Base: show

include("simulations.jl")
include("utils.jl")
include("update.jl")

end

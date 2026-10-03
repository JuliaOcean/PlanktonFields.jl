module Carbon

using KernelAbstractions

using PlanktonKernels.Architectures: device, Architecture
using PlanktonKernels.Grids: AbstractGrid
using PlanktonKernels.Fields: Field, fill_halo_tracer!, init_tracers, initialize_tracer!

using PlanktonFields: PlanktonTracer, PlanktonTracers, CarbonMode, plankton_setup

include("plankton_generation.jl")
include("plankton_growth.jl")
include("plankton_update.jl")

end

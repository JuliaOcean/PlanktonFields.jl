module Carbon

using ..Plankton: plankton_transport!, shape_func_dec

using PlanktonKernels.Architectures: Architecture, device
using PlanktonKernels.Grids: AbstractGrid
using PlanktonKernels.Fields: fill_halo_tracer!, init_tracers, initialize_tracer!

using PlanktonFields: PlanktonTracer

using KernelAbstractions: @kernel, @index

include("plankton_generation.jl")
include("plankton_kernels.jl")
include("plankton_growth.jl")
include("plankton_update.jl")

end

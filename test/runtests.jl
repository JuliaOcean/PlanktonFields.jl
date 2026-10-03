using Test
using PlanktonFields
using PlanktonFields.Simulation: update!
using PlanktonKernels: CPU, RectilinearGrid, interior
using PlanktonKernels.Biogeochemistry: bgc_tracer_default_init, bgc_tracer_names
using PlanktonKernels.Grids: volume
using PlanktonFields.Plankton: phyto_params_default, zoo_params_default, update_phyto_params, update_zoo_params

include("test_helpers.jl")

@testset "PlanktonFields" begin
    @testset "Unit tests" begin
        include("model_test.jl")
        include("parameter_test.jl")
        include("simulation_test.jl")
        include("diagnostics_test.jl")
        include("output_test.jl")
    end
    @testset "Example tests" begin
        include("test_0D_carbon_mode.jl")
        include("test_1D_carbon_mode.jl")
        include("test_2D_carbon_mode.jl")
        include("test_3D_carbon_mode.jl")
    end
end

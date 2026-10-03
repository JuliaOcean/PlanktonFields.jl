module PlanktonFields

if VERSION < v"1.11"
    error("This version of PlanktonFields.jl requires Julia v1.11 or newer.")
end

export
    # Types and structs
    CarbonMode, QuotaMode, plankton_setup,
    
    # Model
    PlanktonModel, model_options,

    # Simulation
    PlanktonSimulation, update!, set_vels_fields!, set_PARF_fields!, set_temp_fields!,

    # Mode defaults and model parameter updates
    phyto_params_default, update_phyto_params,
    zoo_params_default, update_zoo_params,
    default_PARF, default_temp,

    # Output
    PlanktonDiagnostics, PlanktonOutputWriter


using PlanktonKernels.Architectures: Architecture
using PlanktonKernels.Grids: AbstractGrid
using PlanktonKernels.Fields: Field

include("model_structs.jl")
include("Plankton/Plankton.jl")
include("Model/Model.jl")
include("Diagnostics/Diagnostics.jl")
include("Output/Output.jl")
include("Simulation/Simulation.jl")

using .Model
using .Plankton
using .Simulation

end

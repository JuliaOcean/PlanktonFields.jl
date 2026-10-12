module Output

import PlanktonFields: PlanktonOutputWriter

using JLD2: jldopen
using Printf: @sprintf

using PlanktonKernels.Fields: interior

using PlanktonFields: PlanktonModel, PlanktonDiagnostics

import Base: show

include("output_writers.jl")
include("write_outputs.jl")

end

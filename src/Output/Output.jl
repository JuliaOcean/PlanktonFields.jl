module Output

import PlanktonFields: PlanktonOutputWriter

export PlanktonOutputWriter

using JLD2
using Printf: @sprintf

using PlanktonKernels.Fields: interior

using PlanktonFields.Diagnostics
using PlanktonFields.Model
using PlanktonFields: PlanktonModel, PlanktonDiagnostics, CarbonMode, QuotaMode

import Base: show


include("output_writers.jl")
include("write_outputs.jl")

end

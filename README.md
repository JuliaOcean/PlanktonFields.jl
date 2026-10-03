# PlanktonFields.jl

**PlanktonFields.jl** is a Julia framework for Eulerian plankton ecosystem modeling using gridded state variables and continuous differential equations. Plankton populations and associated biogeochemical quantities are represented as spatial fields whose temporal evolution is described through coupled differential equations for processes such as growth, nutrient uptake, mortality, and elemental cycling.

PlanktonFields.jl is built on **PlanktonKernels.jl**, which provides the shared computational infrastructure for plankton ecosystem models, including computational architectures, spatial grids, fields, transport operators, tracers, and related numerical utilities. By relying on the same numerical foundation, PlanktonFields.jl can remain focused on the formulation and solution of Eulerian plankton ecosystem models.

Together with **PlanktonIndividuals.jl**, PlanktonFields.jl forms part of a unified modeling framework for studying plankton ecosystems using complementary representations. PlanktonIndividuals.jl explicitly represents discrete plankton cells and their individual states and trajectories, whereas PlanktonFields.jl represents plankton populations as continuous fields governed by differential equations. Both frameworks share common infrastructure through PlanktonKernels.jl, facilitating consistent model development and direct comparisons between individual-based and Eulerian approaches.

## Quick start

```julia
using PlanktonFields
using PlanktonKernels

grid = RectilinearGrid(size=(8, 8, 4), x=(0, 100), y=(0, 100), z=(0, -20))
model = PlanktonModel(CPU(), grid)
diags = PlanktonDiagnostics(model; tracer=(:DIC, :NO3))
sim = PlanktonSimulation(model; ΔT=60, iterations=10, diags=diags)
```

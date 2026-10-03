function model_options(;kc::Float64 = 0.04, kw::Float64 = 0.046)
    return ModelOpts(kc, kw)
end

"""
    PlanktonModel(arch::Architecture, grid::AbstractGrid;
                  FT = Float32,
                  mode::AbstractMode = CarbonMode(),
                  options::ModelOpts = model_options(),
                  bgc_params = nothing, 
                  tracer_initial = bgc_tracer_default_init(),
                  plankton = nothing,
                  t::Real = 0.0f0,
                  )

Generate a `PlanktonModel` data structure. 

Keyword Arguments (Required)
============================
- `arch` : `CPU()` or `GPU()`. Computer architecture being used to run the model.
- `grid` : a `AbstractGrid` structure. Discrete grid for the model (resolution and geometry).

Keyword Arguments (Optional)
============================
- `options`: Model settings from `model_options()`, stored as `model.options`.
                Includes light attenuation coefficients `kc` and `kw`.
- `FT`: Floating point data type. Default: `Float32`.
- `mode` : Phytoplankton physiology mode, choose among CarbonMode() or CarbonMode().
- `bgc_params` : Parameter set for biogeochemical processes modeled in the model, use default if `nothing`, 
                    use `Dict` to update parameters, the format and names of parameters can be found by running `bgc_params_default()`.
- `tracer_initial` : The source of initial conditions of tracer fields, should be either a `NamedTuple` 
                    or a `Dict` containing the file paths pointing to the files of nutrient initial conditions.
- `plankton`: `plankton_setup(params=..., Nsp=...)` for CarbonMode species; `nothing` creates no plankton.
- `t` : Model time, start from 0 by default, in second.
"""
function PlanktonModel(arch::Architecture, grid::AbstractGrid;
                       FT = Float32,
                       mode::AbstractMode = CarbonMode(),
                       options::ModelOpts = model_options(),
                       bgc_params = nothing, 
                       tracer_initial = bgc_tracer_default_init(),
                       plankton = nothing,
                       t::Real = 0.0f0,
                       )

    isfunctional(arch) || throw(ArgumentError("The requested architecture is unavailable"))

    FT <: AbstractFloat || throw(ArgumentError("FT must be an AbstractFloat type"))

    isfinite(FT(t)) || throw(ArgumentError("Model time must be finite"))

    bgc_params === nothing || bgc_params isa Dict ||
        throw(ArgumentError("Biogeochemical parameters must be nothing or a Dict"))
    bgc_params_final = update_bgc_params(bgc_params === nothing ? Dict() : bgc_params, FT)

    grid_d = replace_grid_storage(arch, grid)

    mode isa CarbonMode || throw(ArgumentError("Only CarbonMode is implemented"))

    arch isa CPU || throw(ArgumentError("CarbonMode currently supports CPU execution"))

    for value in (options.kc, options.kw)
        isfinite(value) && value >= 0 || throw(ArgumentError("Light attenuation must be finite and nonnegative"))
    end

    if plankton === nothing
        phyto_Nsp = 1
        zoo_Nsp = 1
        phyto_init = phyto_tracer_default_init(mode)
        zoo_init = zoo_tracer_default_init(mode)
        phyto_params = phyto_params_default(phyto_Nsp, mode)
        zoo_params = zoo_params_default(zoo_Nsp, mode)
        plankton = plankton_setup(phyto_params, zoo_params, phyto_init, zoo_init, phyto_Nsp, zoo_Nsp)
    elseif plankton isa plankton_setup
        plankton.phyto_params = update_phyto_params(plankton.phyto_params, FT; N=plankton.phyto_Nsp, mode=mode)
        plankton.zoo_params = update_zoo_params(plankton.zoo_params, FT; N=plankton.zoo_Nsp, mode=mode)
    else
        throw(ArgumentError("plankton must be nothing or plankton_setup"))
    end
    plankton_tracers = generate_plankton_tracers(plankton, arch, FT, grid_d, mode)

    plankton_tracer_names = keys(plankton_tracers.phytos)

    bgc_tracers = generate_bgc_tracers(arch, grid_d, tracer_initial, FT)

    ts = timestepper(plankton_tracer_names, arch, FT, grid_d)

    iteration  = 0

    model = PlanktonModel(arch, options, FT, FT(t), iteration, plankton_tracers, bgc_tracers, grid_d, bgc_params_final, ts, mode)

    return model
end

function show(io::IO, model::PlanktonModel)
    Nsp = length(model.plankton.phytos)

    print(io, "PlanktonModel:\n",
              "├── floating point data type: $(model.FT)\n",
              "├── grid: $(short_show(model.grid))\n",
              "├── $(model.mode) is selected for phytoplankton physiology\n",
              "└── phytoplankton: $(Nsp) species\n")
end

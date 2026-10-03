function timestepper(plankton_tracer_names, arch::Architecture, FT::DataType, g::AbstractGrid)
    vel₀ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel½ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))
    vel₁ = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))

    bgc_Gcs = init_tracers(arch, g, bgc_tracer_names, FT)
    bgc_tracer_tmp = init_tracers(arch, g, bgc_tracer_names, FT)
    plk = init_tracers(arch, g, bgc_tracer_names, FT)
    plankton_Gcs = init_tracers(arch, g, plankton_tracer_names, FT)
    plankton_tracer_tmp = init_tracers(arch, g, plankton_tracer_names, FT)

    par = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    par₀= zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    Chl = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    flux_sink = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    temp = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(FT, g.Nx, g.Ny) |> array_type(arch)

    ts = timestepper(bgc_Gcs, bgc_tracer_tmp, plankton_Gcs, plankton_tracer_tmp, vel₀, vel½, vel₁, PARF, temp, flux_sink, plk, par, par₀, Chl)

    return ts
end

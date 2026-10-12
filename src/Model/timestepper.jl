function timestepper(plankton_tracers, arch::Architecture, FT::DataType, g::AbstractGrid, palat)
    vel = (u = Field(arch, g, FT), v = Field(arch, g, FT), w = Field(arch, g, FT))

    bgc_Gcs = init_tracers(arch, g, bgc_tracer_names, FT)
    bgc_tmp = init_tracers(arch, g, bgc_tracer_names, FT)
    plk = init_tracers(arch, g, bgc_tracer_names, FT)

    phyto_tracer_names = keys(plankton_tracers.phytos.sp1.fields)
    phyto_Gcs = init_tracers(arch, g, phyto_tracer_names, FT)
    phyto_tmp = init_tracers(arch, g, phyto_tracer_names, FT)
    
    zoo_tracer_names = keys(plankton_tracers.zoos.sp1.fields)
    zoo_Gcs = init_tracers(arch, g, zoo_tracer_names, FT)
    zoo_tmp = init_tracers(arch, g, zoo_tracer_names, FT)

    par = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    par₀= zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    Chl = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    f_sink = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    temp = zeros(FT, g.Nx+g.Hx*2, g.Ny+g.Hy*2, g.Nz+g.Hz*2) |> array_type(arch)
    PARF = zeros(FT, g.Nx, g.Ny) |> array_type(arch)

    ts = timestepper(bgc_Gcs, bgc_tmp, phyto_Gcs, phyto_tmp, 
                     zoo_Gcs, zoo_tmp, vel, PARF, temp, 
                     f_sink, plk, par, par₀, Chl, palat)

    return ts
end

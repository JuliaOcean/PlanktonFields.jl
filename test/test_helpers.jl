test_grid(size=(2,2,2)) = RectilinearGrid(; size, x=(0,2size[1]), y=(0,2size[2]), z=(0,-2size[3]))
function test_bgc_initial()
    defaults = bgc_tracer_default_init()
    return NamedTuple{keys(defaults)}(Tuple((init=v.init, rand_noise=0.0) for v in defaults))
end
test_model(; size=(2,2,2), kwargs...) = PlanktonModel(CPU(), test_grid(size); tracer_initial=test_bgc_initial(), kwargs...)
function field_mass(field, grid)
    return sum(field.data[i+grid.Hx,j+grid.Hy,k+grid.Hz] *
               volume(i+grid.Hx,j+grid.Hy,k+grid.Hz,grid)
               for i in 1:grid.Nx, j in 1:grid.Ny, k in 1:grid.Nz)
end
function carbon_mass(model)
    mass = sum(field_mass(model.bgc_tracers[n],model.grid) for n in (:DIC,:DOC,:POC))
    for species in (model.plankton.phytos,model.plankton.zoos), plank in species
        mass += field_mass(plank.data.C,model.grid)
    end
    return mass
end
function test_carbon_example(size; moving=false)
    model = test_model(; size, FT=Float64, t=120.0)
    initial_mass = carbon_mass(model)
    nx,ny,nz=size
    # Periodic horizontal flow; bounded top and bottom remain impermeable.
    vels = moving ? (u=fill(1e-4,nx,ny,nz,8), v=zeros(nx,ny,nz,8), w=zeros(nx,ny,nz+1,8)) : (;)
    sim = PlanktonSimulation(model; ΔT=30, iterations=4, vels, diags=test_diagnostics(model),
        PARF=fill(100.0,nx,ny,1),temp=fill(25.0,nx,ny,nz,1))
    update!(sim)
    @test model.iteration == 4
    @test model.t == 240.0
    @test isapprox(carbon_mass(model),initial_mass; rtol=1e-10,atol=1e-10)
    for field in model.bgc_tracers
        @test all(isfinite,interior(field,model.grid))
        @test minimum(interior(field,model.grid)) >= -1e-12
    end
    for species in (model.plankton.phytos,model.plankton.zoos), plank in species
        @test all(isfinite,interior(plank.data.C,model.grid))
        @test minimum(interior(plank.data.C,model.grid)) >= 0
    end
end

function test_diagnostics(model; kwargs...)
    return PlanktonDiagnostics(model;
        phytoplankton=keys(model.plankton.phytos.sp1.data),
        zooplankton=keys(model.plankton.zoos.sp1.data), kwargs...)
end

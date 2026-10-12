@testset "Model options" begin
    options=model_options(kc=0.08,kw=0.02)
    @test options.kc == 0.08
    @test options.kw == 0.02
end
@testset "Model initialization ($FT)" for FT in (Float32,Float64)
    model=test_model(; FT,bgc_params=Dict("κh"=>0.2),t=5)
    @test keys(model.bgc_tracers) == bgc_tracer_names
    @test all(==(FT(0.8)),interior(model.bgc_tracers.NO3,model.grid))
    @test model.bgc_params["κh"] === FT(0.2)
    @test model.t === FT(5)
    @test model.iteration == 0
    @test keys(model.plankton.phytos) == (:sp1,)
    @test keys(model.plankton.zoos) == (:sp1,)
    @test eltype(model.plankton.phytos.sp1.fields.C.data) === FT
    @test eltype(model.plankton.zoos.sp1.fields.C.data) === FT
    @test model.timestepper.bgc_tmp.NO3.data !== model.bgc_tracers.NO3.data
    @test model.timestepper.plk.NO3.data !== model.timestepper.bgc_Gcs.NO3.data
    @test !isempty(sprint(show,model))
end
@testset "Invalid model configuration" begin
    @test_throws ArgumentError test_model(t=NaN)
    @test_throws ArgumentError test_model(FT=Int)
    @test_throws ArgumentError test_model(bgc_params=42)
    @test_throws ArgumentError test_model(mode=QuotaMode())
    @test_throws ArgumentError test_model(options=model_options(kc=-1.0))
end

@testset "Plankton transport increments" begin
    model=test_model(size=(4,2,2),FT=Float64)
    grid=model.grid
    dt=3.0
    model.timestepper.vel.u.data .= 0.01
    for (plank,work,tmp) in
        ((model.plankton.phytos.sp1,model.timestepper.phyto_Gcs,model.timestepper.phyto_tmp),
         (model.plankton.zoos.sp1,model.timestepper.zoo_Gcs,model.timestepper.zoo_tmp))
        plank.fields.C.data .= 1
        plank.fields.C.data[grid.Hx+2,grid.Hy+1,grid.Hz+1]=2
        PlanktonKernels.Fields.fill_halo_tracer!(plank.fields,grid)
        before=copy(plank.fields.C.data)
        PlanktonKernels.Transport.tracer_advection!(plank.fields,tmp,work,model.timestepper.vel,grid,dt,CPU())
        PlanktonKernels.Transport.tracer_diffusion!(work,CPU(),grid,plank.fields,plank.p.κh,plank.p.κh,plank.p.κv,dt)
        expected=before .+ work.C.data
        @test any(!iszero,interior(work.C,grid))
        PlanktonFields.Plankton.plankton_transport!(model.plankton,model.timestepper,grid,dt,CPU())
        @test interior(plank.fields.C,grid) ≈ interior(expected,grid)
    end
end

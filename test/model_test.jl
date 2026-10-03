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
    @test eltype(model.plankton.phytos.sp1.data.C.data) === FT
    @test eltype(model.plankton.zoos.sp1.data.C.data) === FT
    @test model.timestepper.bgc_tracer_tmp.NO3.data !== model.bgc_tracers.NO3.data
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

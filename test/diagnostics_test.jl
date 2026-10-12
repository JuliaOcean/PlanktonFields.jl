@testset "Diagnostics" begin
    model=test_model()
    diags=test_diagnostics(model;tracer=(:DIC,:PIFe,:POFe),iteration_interval=2)
    @test diags.iteration_interval == 2
    @test all(n -> haskey(diags.bgc_tracers,n),(:DIC,:PIFe,:POFe,:PAR))
    @test keys(diags.phytos) == keys(model.plankton.phytos)
    @test keys(diags.zoos) == keys(model.plankton.zoos)
    @test all(iszero,diags.bgc_tracers.DIC)
    @test !isempty(sprint(show,diags))
end

@testset "Default diagnostics" begin
    model=test_model()
    diags=PlanktonDiagnostics(model)
    @test keys(diags.phytos) == (:sp1,)
    @test keys(diags.zoos) == (:sp1,)
end

@testset "State and process diagnostics through timestepping" begin
    model=test_model()
    phyto=model.plankton.phytos.sp1
    zoo=model.plankton.zoos.sp1
    phyto.fields.C.data .= 3
    phyto.fields.Chl.data .= 0.3
    zoo.fields.C.data .= 5
    diags=PlanktonDiagnostics(model;phytoplankton=(:C,:graz),zooplankton=(:C,:G))
    expected_phyto=zeros(size(phyto.fields.C.data))
    expected_zoo=zeros(size(zoo.fields.C.data))
    expected_graz=zeros(size(phyto.procs.graz.data))
    expected_G=zeros(size(zoo.procs.G.data))
    for _ in 1:2
        PlanktonFields.Model.TimeStep!(model,1.0f0,diags)
        expected_phyto .+= phyto.fields.C.data
        expected_zoo .+= zoo.fields.C.data
        expected_graz .+= phyto.procs.graz.data
        expected_G .+= zoo.procs.G.data
    end
    @test diags.phytos.sp1.C ≈ expected_phyto
    @test diags.zoos.sp1.C ≈ expected_zoo
    @test diags.phytos.sp1.graz ≈ expected_graz
    @test diags.zoos.sp1.G ≈ expected_G
    @test any(>(0),expected_graz)
end

@testset "State and process diagnostic output" begin
    model=test_model()
    diags=PlanktonDiagnostics(model;phytoplankton=(:C,:PS),zooplankton=(:G,))
    # Supply accumulated samples directly to isolate output from physiology.
    diags.phytos.sp1.C .= 6
    diags.phytos.sp1.PS .= 14
    diags.zoos.sp1.G .= 22
    @test keys(diags.zoos.sp1) == (:G,)
    @test_throws ArgumentError PlanktonDiagnostics(model;phytoplankton=(:unknown,))
    mktempdir() do dir
        writer=PlanktonOutputWriter(;dir,save_diags=true)
        model.iteration=2
        diags.iteration_interval=2
        PlanktonFields.Output.write_output!(writer,model,diags)
        PlanktonFields.Output.jldopen(writer.diags_file,"r") do file
            @test all(==(3),file["timeseries/phytoplankton/sp1/C/2"])
            @test all(==(7),file["timeseries/phytoplankton/sp1/PS/2"])
            @test all(==(11),file["timeseries/zooplankton/sp1/G/2"])
        end
        @test all(iszero,diags.phytos.sp1.PS)
        @test all(iszero,diags.phytos.sp1.C)
        @test all(iszero,diags.zoos.sp1.G)
    end
end

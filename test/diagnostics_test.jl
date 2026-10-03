@testset "Diagnostics" begin
    model=test_model()
    diags=test_diagnostics(model;tracer=(:DIC,:PIFe,:POFe),iteration_interval=2)
    @test diags.iteration_interval == 2
    @test all(n -> haskey(diags.bgc_tracers,n),(:DIC,:PIFe,:POFe,:PAR))
    @test keys(diags.phyto_tracers) == keys(model.plankton.phytos)
    @test keys(diags.zoo_tracers) == keys(model.plankton.zoos)
    @test all(iszero,diags.bgc_tracers.DIC)
    @test !isempty(sprint(show,diags))
end

@testset "Default diagnostics" begin
    model=test_model()
    diags=PlanktonDiagnostics(model)
    @test keys(diags.phyto_tracers) == (:sp1,)
    @test keys(diags.zoo_tracers) == (:sp1,)
end

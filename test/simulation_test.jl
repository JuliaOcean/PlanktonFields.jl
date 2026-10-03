@testset "Simulation and diagnostic configuration" begin
    grid = RectilinearGrid(size=(2, 2, 2), x=(0, 2), y=(0, 2), z=(0, -2))
    model = test_model()
    diags = test_diagnostics(model; tracer=(:DIC, :NO3), iteration_interval=2)
    sim = PlanktonSimulation(model; ΔT=60, iterations=10, diags=diags)
    @test sim.model === model
    @test sim.ΔT === 60.0f0
    @test sim.iterations == 10
    @test sim.diags === diags
    @test sim.output_writer === nothing
    @test_throws ArgumentError test_diagnostics(model; tracer=(:unknown,))
    @test_throws ArgumentError test_diagnostics(model; tracer=(:DIC, :DIC))
    @test_throws ArgumentError test_diagnostics(model; iteration_interval=0)
    @test_throws ArgumentError PlanktonSimulation(model; ΔT=0, iterations=1)
    @test_throws ArgumentError PlanktonSimulation(model; ΔT=Inf, iterations=1)
    @test_throws ArgumentError PlanktonSimulation(model; ΔT=1, iterations=-1)
end

@testset "Forcing coverage and precision" begin
    model=test_model()
    diags=test_diagnostics(model)
    sim=PlanktonSimulation(model; ΔT=60,iterations=1,diags)
    @test eltype(sim.input.temp) === Float32
    @test eltype(sim.input.PARF) === Float32
    @test !isempty(sprint(show,sim))
    @test_throws ArgumentError PlanktonSimulation(model; ΔT=60,iterations=61,diags,
        PARF=zeros(Float32,2,2,1))
    @test_throws ArgumentError PlanktonSimulation(model; ΔT=60,iterations=1,diags,ΔT_temp=0)
end

using JLD2: jldopen
using PlanktonFields.Output: write_output!

@testset "Output writer configuration" begin
    mktempdir() do root
        writer=PlanktonOutputWriter(dir=joinpath(root,"results"),save_diags=true)
        @test isdir(writer.filepath)
        @test endswith(writer.diags_file,"diags.jld2")
        @test !isempty(sprint(show,writer))
        @test_throws ArgumentError PlanktonOutputWriter(dir=joinpath(root,"invalid"),max_filesize=0)
    end
end
@testset "Output averaging and file rotation" begin
    model=test_model()
    diags=test_diagnostics(model;tracer=(:NO3,),iteration_interval=2)
    mktempdir() do dir
        writer = PlanktonOutputWriter(; dir, save_diags=true, max_filesize=1)
        @test !isempty(sprint(show, writer))
        model.iteration = 2
        diags.bgc_tracers.NO3 .= 4
        write_output!(writer, model, diags)
        @test isfile(writer.diags_file)
        jldopen(writer.diags_file, "r") do file
            @test all(file["timeseries/NO3/2"] .== 2)
        end
        @test all(iszero, diags.bgc_tracers.NO3)
        model.iteration = 4
        diags.bgc_tracers.NO3 .= 6
        write_output!(writer, model, diags)
        @test writer.part_diags == 2
        @test isfile(joinpath(dir, "diags_part1.jld2"))
        jldopen(writer.diags_file, "r") do file
            @test all(file["timeseries/NO3/4"] .== 3)
        end
    end
end

function palatability_setup(palat; nphyto=1, nzoo=1)
    mode=CarbonMode()
    phyto_names=ntuple(i -> Symbol("sp",i),nphyto)
    zoo_names=ntuple(i -> Symbol("sp",i),nzoo)
    phyto_init=NamedTuple{phyto_names}(ntuple(_ -> (Chl=(init=0.1,rand_noise=0.0),C=(init=1.0,rand_noise=0.0)),nphyto))
    zoo_init=NamedTuple{zoo_names}(ntuple(_ -> (C=(init=1.0,rand_noise=0.0),),nzoo))
    return plankton_setup(phyto_params_default(nphyto,mode),zoo_params_default(nzoo,mode),
                         phyto_init,zoo_init,nphyto,nzoo,palat)
end

@testset "Palatability" begin
    @testset "Precision and independent storage ($FT)" for FT in (Float32,Float64)
        source=[0.0 0.25 1.0; 1.0 0.75 0.0]
        setup=palatability_setup(source;nphyto=2,nzoo=3)
        model=test_model(;FT,plankton=setup)
        @test model.timestepper.palat isa Matrix{FT}
        @test size(model.timestepper.palat) == (2,3)
        @test model.timestepper.palat == FT.(source)
        @test setup.palat === source
        @test model.timestepper.palat !== source
        source[1,1]=0.5
        @test model.timestepper.palat[1,1] == 0
        default=test_model(;FT)
        @test default.timestepper.palat == ones(FT,1,1)
        @test eltype(default.timestepper.palat) === FT
    end
    @testset "Invalid entries and dimensions" begin
        for value in (-0.1,1.1,NaN,Inf,-Inf,nextfloat(1.0),-eps(Float64),1+0im,"0.5")
            @test_throws ArgumentError test_model(plankton=palatability_setup(fill(value,1,1)))
        end
        for array in (ones(2,1),ones(1,2),ones(1),ones(1,1,1))
            @test_throws ArgumentError test_model(plankton=palatability_setup(array))
        end
    end
    @testset "Array input conversion" begin
        for array in (ones(Int,1,1),reshape([true],1,1),view(fill(0.5,2,2),1:1,1:1))
            model=test_model(plankton=palatability_setup(array))
            @test model.timestepper.palat isa Matrix{Float32}
            @test model.timestepper.palat == Float32.(array)
        end
    end
end

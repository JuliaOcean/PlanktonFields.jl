@testset "Parameter defaults and overrides" begin
    for defaults in (phyto_params_default,zoo_params_default), n in (1,3)
        params=defaults(n,CarbonMode())
        @test all(v -> length(v)==n,values(params))
        @test all(v -> all(isfinite,v),values(params))
        firstkey=first(keys(params))
        original=params[firstkey][1]
        params[firstkey][1]=-1
        @test defaults(n,CarbonMode())[firstkey][1] == original
    end
    for (defaults,update) in ((phyto_params_default,update_phyto_params),(zoo_params_default,update_zoo_params)), FT in (Float32,Float64)
        override=Dict("k_RS"=>[2e-6,3e-6])
        saved=deepcopy(override)
        params=update(override,FT;N=2,mode=CarbonMode())
        @test params["k_RS"] == FT[2e-6,3e-6]
        @test override == saved
        @test keys(params) == keys(defaults(2,CarbonMode()))
        @test_throws ArgumentError update(Dict("unknown"=>[1.0]),FT;mode=CarbonMode())
    end
end
@testset "Default parameters construct each plankton type" begin
    grid=test_grid()
    carbon=PlanktonFields.Plankton.Carbon
    @testset "$kind" for (kind,defaults,construct) in
        ((:phyto,phyto_params_default,carbon.construct_phytoplankton),
         (:zoo,zoo_params_default,carbon.construct_zooplankton))
        plank=construct(CPU(),1,defaults(1,CarbonMode()),Float64,grid)
        @test eltype(plank.fields.C.data) === Float64
        @test all(iszero,plank.fields.C.data)
        @test keys(plank.fields) == (kind == :phyto ? (:Chl,:C) : (:C,))
        @test keys(plank.procs) == (kind == :phyto ? (:PS,:VNO3,:VNH4,:VPO4,:BS,:SChl,:RS,:mort,:graz) : (:G,:RS,:mort))
        @test all(f -> eltype(f.data) === Float64 && all(iszero,f.data), plank.procs)
        @test plank.fields.C.data !== plank.procs.RS.data
    end
end

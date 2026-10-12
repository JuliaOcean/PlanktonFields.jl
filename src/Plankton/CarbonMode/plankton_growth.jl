function phyto_growth!(phytos, bgc_tracers, timestepper, grid, ΔT, arch)
    par = timestepper.par
    Temp = timestepper.temp
    for sp in eachindex(phytos)
        phyto = phytos[sp]
        p = phytos[sp].p
        calc_PS!(phyto, Temp, par, p, grid, arch)
        calc_VN!(phyto, bgc_tracers, Temp, p, ΔT, grid, arch)
        calc_BS!(phyto, p, grid, arch)
        calc_SChl!(phyto, par, Temp, p, grid, arch)
        calc_respir!(phyto, Temp, p, grid, arch)
        calc_mort!(phyto, p, grid, arch)
    end
    return nothing
end

function zoo_growth!(zoos, timestepper, grid, arch)
    Temp = timestepper.temp
    for sp in eachindex(zoos)
        zoo = zoos[sp]
        p = zoos[sp].p
        calc_respir!(zoo, Temp, p, grid, arch)
        calc_quadratic_mort!(zoo, p, grid, arch)
    end
    return nothing
end

function grazing!(plankton, timestepper, ΔT)
    palat = timestepper.palat
    plk = timestepper.plk
    graz = timestepper.f_sink
    for i in axes(palat, 1), j in axes(palat, 2)
        phyto = plankton.phytos[i]
        zoo = plankton.zoos[j]

        ##### calculate grazing
        @inbounds graz .= phyto.fields.C.data .* zoo.fields.C.data .* palat[i,j] .* zoo.p.Gmax
        @inbounds phyto.procs.graz.data .+= graz
        @inbounds zoo.procs.G.data .+= graz .* zoo.p.Gfrac

        ##### grazing release: sloopy feeding
        @inbounds plk.DOC.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* phyto.p.grazFracC .* ΔT
        @inbounds plk.POC.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* 
                                        (1.0f0 - phyto.p.grazFracC) .* ΔT

        @inbounds plk.DON.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* phyto.p.R_NC .* 
                                        phyto.p.grazFracC .* ΔT
        @inbounds plk.PON.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* phyto.p.R_NC .* 
                                        (1.0f0 - phyto.p.grazFracC) .* ΔT

        @inbounds plk.DOP.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* phyto.p.R_PC .* 
                                        phyto.p.grazFracC .* ΔT
        @inbounds plk.POP.data .+= graz .* (1.0f0 - zoo.p.Gfrac) .* phyto.p.R_PC .* 
                                        (1.0f0 - phyto.p.grazFracC) .* ΔT
    end
    return nothing
end

function reset_grazing!(plankton)
    phytos = plankton.phytos
    zoos = plankton.zoos
    for sp in eachindex(phytos)
        phytos[sp].procs.graz.data .= 0.0f0
    end
    for sp in eachindex(zoos)
        zoos[sp].procs.G.data .= 0.0f0
    end
    return nothing
end

function apply_phytos_tendencies!(phytos, phyto_tmp, ΔT, grid::AbstractGrid)
    for sp in eachindex(phytos)
        phyto = phytos[sp]
        @inbounds phyto_tmp.C.data .= phyto.fields.C.data
        @inbounds phyto_tmp.Chl.data .= phyto.fields.Chl.data

        @inbounds phyto.fields.C.data .+= phyto.procs.BS.data .* ΔT
        @inbounds phyto.fields.C.data .-= phyto.procs.RS.data .* ΔT
        @inbounds phyto.fields.C.data .-= phyto.procs.mort.data .* ΔT
        @inbounds phyto.fields.C.data .-= phyto.procs.graz.data .* ΔT

        @inbounds phyto.fields.Chl.data .+= phyto.procs.SChl.data .* ΔT
        @inbounds phyto.fields.Chl.data .-= (phyto.procs.mort.data .+ 
                                             phyto.procs.graz.data) .*
                                             phyto_tmp.Chl.data ./ 
                                             max.(1.0f-30, phyto_tmp.C.data) .* ΔT
        fill_halo_tracer!(phyto.fields, grid)
    end
    return nothing
end

function apply_zoos_tendencies!(zoos, ΔT, grid::AbstractGrid)
    for sp in eachindex(zoos)
        zoo = zoos[sp]
        @inbounds zoo.fields.C.data .+= zoo.procs.G.data .* ΔT
        @inbounds zoo.fields.C.data .-= zoo.procs.RS.data .* ΔT
        @inbounds zoo.fields.C.data .-= zoo.procs.mort.data .* ΔT
        fill_halo_tracer!(zoo.fields, grid)
    end
    return nothing
end

function calc_consume_loss!(plankton, plk, ΔT)
    for sp in eachindex(plankton.phytos)
        phyto = plankton.phytos[sp]
        p = plankton.phytos[sp].p
        plk.DIC.data .+= (phyto.procs.RS.data .- phyto.procs.BS.data) .* ΔT
        plk.NH4.data .-= phyto.procs.BS.data .* p.R_NC .* ΔT .* phyto.procs.VNH4.data ./ 
                            max.(1.0f-30, phyto.procs.VNH4.data .+ phyto.procs.VNO3.data)
        plk.NO3.data .-= phyto.procs.BS.data .* p.R_NC .* ΔT.* phyto.procs.VNO3.data ./ 
                            max.(1.0f-30, phyto.procs.VNH4.data .+ phyto.procs.VNO3.data)
        plk.PO4.data .-= phyto.procs.BS.data .* p.R_PC .* ΔT
        plk.DOC.data .+= phyto.procs.mort.data .* p.mortFracC .* ΔT
        plk.POC.data .+= phyto.procs.mort.data .* (1.0f0 .- p.mortFracC) .* ΔT
        plk.DON.data .+= phyto.procs.mort.data .* p.R_NC .* p.mortFracC .* ΔT
        plk.PON.data .+= phyto.procs.mort.data .* p.R_NC .* (1.0f0 .- p.mortFracC) .* ΔT
        plk.DOP.data .+= phyto.procs.mort.data .* p.R_PC .* p.mortFracC .* ΔT
        plk.POP.data .+= phyto.procs.mort.data .* p.R_PC .* (1.0f0 .- p.mortFracC) .* ΔT
    end

    for sp in eachindex(plankton.zoos)
        zoo = plankton.zoos[sp]
        p = plankton.zoos[sp].p
        plk.DIC.data .+= zoo.procs.RS.data .* ΔT
        plk.DOC.data .+= zoo.procs.mort.data .* p.mortFracC .* ΔT
        plk.POC.data .+= zoo.procs.mort.data .* (1.0f0 .- p.mortFracC) .* ΔT
        plk.DON.data .+= zoo.procs.mort.data .* p.R_NC .* p.mortFracC .* ΔT
        plk.PON.data .+= zoo.procs.mort.data .* p.R_NC .* (1.0f0 .- p.mortFracC) .* ΔT
        plk.DOP.data .+= zoo.procs.mort.data .* p.R_PC .* p.mortFracC .* ΔT
        plk.POP.data .+= zoo.procs.mort.data .* p.R_PC .* (1.0f0 .- p.mortFracC) .* ΔT
    end
    return nothing
end
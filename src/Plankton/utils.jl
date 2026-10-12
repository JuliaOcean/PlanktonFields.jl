##### calculate PAR field based on Chla field and depth
@kernel function calc_par_kernel!(par, Chl, PARF, g::AbstractGrid, kc, kw, ki)
    i, j = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = ki + g.Hz
    atten = (Chl[ii,jj,kk]/volume(ii, jj, kk, g) * kc + kw) * ΔzF(ii, jj, kk, g)
    par[ii,jj,kk] = PARF[i,j] * (1.0f0 - exp(-atten)) / atten
    PARF[i,j] = PARF[i,j] * exp(-atten)
end
function calc_par!(par, arch::Architecture, Chl, PARF, g::AbstractGrid, kc, kw, ki)
    kernel! = calc_par_kernel!(device(arch), (16,16), (g.Nx, g.Ny))
    kernel!(par, Chl, PARF, g, kc, kw, ki)
    return nothing
end

##### calculate the advection and diffusiton of plankton tracers
function plankton_transport!(plankton::PlanktonTracers, timestepper, grid::AbstractGrid, ΔT, arch::Architecture)
    phyto_tmp = timestepper.phyto_tmp
    phyto_Gcs = timestepper.phyto_Gcs
    zoo_tmp = timestepper.zoo_tmp
    zoo_Gcs = timestepper.zoo_Gcs
    vel = timestepper.vel

    for sp in plankton.phytos
        zero_fields!(phyto_tmp)
        zero_fields!(phyto_Gcs)
        tracer_advection!(sp.fields, phyto_tmp, phyto_Gcs, vel, grid, ΔT, arch)
        tracer_diffusion!(phyto_Gcs, arch, grid, sp.fields, sp.p.κh, sp.p.κh, sp.p.κv, ΔT)
        ##### apply tendencies to update tracers
        for tr in keys(sp.fields)
            @inbounds sp.fields[tr].data .+= phyto_Gcs[tr].data
        end
        fill_halo_tracer!(sp.fields, grid)
    end

    for sp in plankton.zoos
        zero_fields!(zoo_tmp)
        zero_fields!(zoo_Gcs)
        tracer_advection!(sp.fields, zoo_tmp, zoo_Gcs, vel, grid, ΔT, arch)
        tracer_diffusion!(zoo_Gcs, arch, grid, sp.fields, sp.p.κh, sp.p.κh, sp.p.κv, ΔT)
        ##### apply tendencies to update tracers
        for tr in keys(sp.fields)
            @inbounds sp.fields[tr].data .+= zoo_Gcs[tr].data
        end
        fill_halo_tracer!(sp.fields, grid)
    end
end

##### shape function - decrease from 1.0 to 0.0 while x increase from 0.0 to 1.0
##### sharp decrease near x/xmax = 1.0
@inline function shape_func_dec(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, 1.0f0 - x / xmax))
    reg = fx^pow / (k + fx^pow)
    return reg
end

##### shape function - increase from 0.0 to 1.0 while x increase from 0.0 to 1.0
##### sharp increase near x/xmax = 1.0
@inline function shape_func_inc(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, 1.0f0 - x / xmax))
    reg = fx^pow / (k + fx^pow)
    return 1.0f0 - reg
end

##### shape function - decrease from 1.0 to 0.0 while x increase from 0.0 to 1.0
##### sharp decrease near x/xmax = 0.0
@inline function shape_func_dec_alt(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, x / xmax))
    reg = fx^pow / (k + fx^pow)
    return 1.0f0 - reg
end

##### shape function - increase from 0.0 to 1.0 while x increase from 0.0 to 1.0
##### sharp increase near x/xmax = 0.0
@inline function shape_func_inc_alt(x, xmax, k; pow = 4.0f0)
    fx = max(0.0f0, min(1.0f0, x / xmax))
    reg = fx^pow / (k + fx^pow)
    return reg
end

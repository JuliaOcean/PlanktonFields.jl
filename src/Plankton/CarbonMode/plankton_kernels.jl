##### temperature function for photosynthesis
@inline function tempFunc_PS(T, p)
    x = T - p.Topt; xmax = p.Tmax - p.Topt
    regT = shape_func_dec(x, xmax, 1.0f-1)
    Ea = p.Ea_ref * (p.Topt_ref / max(1.0f-30, p.Topt))
    Δtemp = 1.0f0/(T+273.15f0) - 1.0f0/(p.Topt+273.15f0)
    k = exp(-Ea/8.3145f0*Δtemp) * regT
    return min(1.0f0, k)
end

##### temperature function for metabolic rates
@inline function tempFunc(T, p)
    Ea = p.Ea_ref * (p.Topt_ref / p.Topt)
    Δtemp = 1.0f0/(T+273.15f0) - 1.0f0/(p.Topt+273.15f0)
    k = exp(-Ea/8.3145f0*Δtemp)
    return min(1.0f0, k)
end

##### calculate photosynthesis (mmolC/m³/second)
@kernel function calc_PS_kernel!(PS, C, Chl, Temp, par, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds αI = par[ii, jj, kk] * p.α * p.Φ
    @inbounds PCm = p.PCmax * tempFunc_PS(Temp[ii, jj, kk], p)
    @inbounds PS[ii, jj, kk] = PCm * (1.0f0 - exp(-αI / max(1.0f-30, PCm) * 
                                Chl[ii, jj, kk] / max(1.0f-30, C[ii, jj, kk]))) *
                                C[ii, jj, kk]
end
function calc_PS!(phyto, Temp, par, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_PS_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.PS.data, phyto.fields.C.data, phyto.fields.Chl.data, 
            Temp, par, p, g)
    return nothing
end

##### calculate nutrient uptake (mmolN/m³/second)
@kernel function calc_VN_kernel!(VN, C, N, Temp, VNmax, KsatN, p, ΔT, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds regN = N[ii, jj, kk] / max(1.0f-30, N[ii, jj, kk] + KsatN)
    @inbounds maxVN = N[ii, jj, kk]/ΔT
    @inbounds VN[ii, jj, kk] = VNmax * regN * tempFunc(Temp[ii, jj, kk], p) * C[ii, jj, kk]
    @inbounds VN[ii, jj, kk] = min(VN[ii, jj, kk], maxVN)
end
function calc_VN!(phyto, bgc, Temp, p, ΔT, g::AbstractGrid, arch::Architecture)
    kernel! = calc_VN_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    ##### calculate NO3 uptake
    kernel!(phyto.procs.VNO3.data, phyto.fields.C.data, bgc.NO3.data, 
            Temp, p.VNO3max, p.KsatNO3, p, ΔT, g)
    ##### calculate NH4 uptake
    kernel!(phyto.procs.VNH4.data, phyto.fields.C.data, bgc.NH4.data, 
            Temp, p.VNH4max, p.KsatNH4, p, ΔT, g)
    ##### calculate PO4 uptake
    kernel!(phyto.procs.VPO4.data, phyto.fields.C.data, bgc.PO4.data, 
            Temp, p.VPO4max, p.KsatPO4, p, ΔT, g)
end

##### calculate biosynthesis rate (mmolC/m³/second)
@kernel function calc_BS_kernel!(BS, PS, VNO3, VNH4, VPO4, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds VN = (VNO3[ii, jj, kk] + VNH4[ii, jj, kk]) / p.R_NC
    @inbounds VP = VPO4[ii, jj, kk] / p.R_PC
    @inbounds BS[ii, jj, kk] = min(PS[ii, jj, kk], VN, VP)
end
function calc_BS!(phyto, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_BS_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.BS.data, phyto.procs.PS.data, phyto.procs.VNO3.data,
            phyto.procs.VNH4.data, phyto.procs.VPO4.data, p, g)
end

@kernel function calc_SChl_kernel!(SChl, PS, C, Chl, par, Temp, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds αI = par[ii, jj, kk] * p.α * p.Φ
    @inbounds Chl2C = Chl[ii, jj, kk] / max(1.0f-30, C[ii, jj, kk])
    @inbounds PCm = p.PCmax * tempFunc_PS(Temp[ii, jj, kk], p)
    @inbounds Chl2C_acclim = p.Chl2Cmax / (1.0f0 + p.Chl2Cmax * αI / 
                                            (max(1.0f-30, 2.0f0 * PCm)))
    @inbounds Chl2C_acclim = clamp(Chl2C_acclim, p.Chl2Cmin, p.Chl2Cmax)
    @inbounds SChl[ii, jj, kk] = Chl2C_acclim * PS[ii, jj, kk] * C[ii, jj, kk] +
                                    p.τChl * (Chl2C_acclim - Chl2C) * C[ii, jj, kk]
end
function calc_SChl!(phyto, par, Temp, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_SChl_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.SChl.data, phyto.procs.PS.data, phyto.fields.C.data,
            phyto.fields.Chl.data, par, Temp, p, g)
end

##### calculate respiration rate (mmolC/m³/second)
@kernel function calc_respir_kernel!(RS, C, Temp, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds RS[ii, jj, kk] = C[ii, jj, kk] * p.k_RS * tempFunc(Temp[ii, jj, kk], p)
end
function calc_respir!(phyto, Temp, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_respir_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.RS.data, phyto.fields.C.data, Temp, p, g)
end

##### calculate mortality rate (mmolC/m³/second)
@kernel function calc_mort_kernel!(mort, C, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds mort[ii, jj, kk] = C[ii, jj, kk] * p.k_mort
end
function calc_mort!(phyto, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_mort_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.mort.data, phyto.fields.C.data, p, g)
end

##### calculate quadratic mortality rate (mmolC/m³/second)
@kernel function calc_quadratic_mort_kernel!(mort, C, p, g)
    i, j, k = @index(Global, NTuple)
    ii = i + g.Hx
    jj = j + g.Hy
    kk = k + g.Hz

    @inbounds mort[ii, jj, kk] = C[ii, jj, kk] * C[ii, jj, kk] * p.k_mort
end
function calc_quadratic_mort!(phyto, p, g::AbstractGrid, arch::Architecture)
    kernel! = calc_quadratic_mort_kernel!(device(arch), (16,16), (g.Nx, g.Ny, g.Nz))
    kernel!(phyto.procs.mort.data, phyto.fields.C.data, p, g)
end
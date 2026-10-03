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

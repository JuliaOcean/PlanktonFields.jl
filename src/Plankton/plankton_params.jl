# Individual parameter defaults, dispatched by physiology mode.
function phyto_params_default end
function zoo_params_default end

function generate_n_species_params(N, params)
    p = []
    for key in keys(params)
        pa = (key, fill(params[key][1],N))
        push!(p, pa)
    end
    return Dict(p)
end

"""
    phyto_params_default(N::Int64, mode::AbstractMode)
Generate default phytoplankton parameter values based on `AbstractMode` and species number `N`.
"""
function phyto_params_default(N::Int64, mode::CarbonMode)
    params=Dict(
        "Cquota"    => [1.8e-11], # C quota of phyto cells at size = 1.0
        "Chl2C"     => [0.10],    # Chla:C ratio in phytoplankton (mgChl/mmolC)
        "α"         => [2.0e-2],  # Irradiance absorption coeff (m²/mgChl)
        "Φ"         => [4.0e-5],  # Maximum quantum yield (mmolC/μmol photon)
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea_ref"    => [5.3e4],   # Free energy
        "Topt_ref"  => [26.0],    # reference optimal growth temperature(C)
        "PCmax"     => [4.2e-5],  # Maximum primary production rate (per second)
        "respir"    => [1.2e-6],  # Respiration rate(per second)
        "k_mort"    => [5e-5],    # Probability of cell natural death per second
        "grazFracC" => [0.7],     # Fraction goes into dissolved organic pool
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    zoo_params_default(N::Int64, mode::AbstractMode)
Generate default zooplankton parameter values based on `AbstractMode` and species number `N`.
"""
function zoo_params_default(N::Int64, mode::CarbonMode)
    params=Dict(
        "Cquota"    => [1.8e-11], # C quota of zoo cells at size = 1.0
        "Topt"      => [27.0],    # Optimal temperature for growth (C)
        "Tmax"      => [30.0],    # Maximal temperature for growth (C)
        "Ea_ref"    => [5.3e4],   # Free energy
        "Topt_ref"  => [26.0],    # reference optimal growth temperature(C)
        "Gmax"      => [1.0e-5],  # Maximum grazing rate (per second)
        "respir"    => [1.2e-6],  # Respiration rate(per second)
        "k_mort"    => [5e-5],    # Probability of cell natural death per second
        "mortFracC" => [0.5],     # Fraction goes into dissolved organic pool
    )

    if N == 1
        return params
    else
        return generate_n_species_params(N, params)
    end
end

"""
    update_phyto_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of phytoplankton physiology resolved in the model
"""
function update_phyto_params(tmp::Dict, FT::DataType; N::Int = 1, mode::AbstractMode = CarbonMode())
    parameters = phyto_params_default(N,mode)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: phyto parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end

"""
    update_zoo_params(tmp::Dict, FT::DataType; N::Int64, mode::AbstractMode)
Update parameter values based on a `Dict` provided by user
Keyword Arguments
=================
- `tmp` is a `Dict` containing the parameters needed to be upadated
- `FT`: Floating point data type. Default: `Float32`.
- `N` is a `Int64` indicating the number of species
- `mode` is the mode of zooplankton physiology resolved in the model
"""
function update_zoo_params(tmp::Dict, FT::DataType; N::Int = 1, mode::AbstractMode = QuotaMode())
    parameters = zoo_params_default(N,mode)
    tmp_keys = collect(keys(tmp))
    pkeys = collect(keys(parameters))
    for key in tmp_keys
        if length(findall(x->x==key, pkeys))==0
            throw(ArgumentError("PARAM: zoo parameter not found $key"))
        else
            parameters[key] = FT.(tmp[key])
        end
    end
    return parameters
end
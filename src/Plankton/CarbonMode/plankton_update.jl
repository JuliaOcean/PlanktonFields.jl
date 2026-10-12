function plankton_update!(plankton, bgc_tracers, timestepper, grid, ΔT, arch::Architecture)
    ##### reset accumulated grazing
    reset_grazing!(plankton)

    ##### advection and diffution
    plankton_transport!(plankton, timestepper, grid, ΔT, arch)

    ##### phytoplankton growth
    phyto_growth!(plankton.phytos, bgc_tracers, timestepper, grid, ΔT, arch)

    ##### grazing
    grazing!(plankton, timestepper, ΔT)

    ##### zooplankton growth
    zoo_growth!(plankton.zoos, timestepper, grid, arch)

    ##### apply tendencies
    apply_phytos_tendencies!(plankton.phytos, timestepper.phyto_tmp, ΔT, grid)
    apply_zoos_tendencies!(plankton.zoos, ΔT, grid)

    ##### calclulate bgc tracers interactions
    calc_consume_loss!(plankton, timestepper.plk, ΔT)

    return nothing
end
"""
    burnin(conf::Union{RotationConfig, PatchConfig}, Nt::Int; show_E=false)

Initialise and run simulation for configuration `conf` with `Nt` timesteps. 
Finds equilibrium position to then be used for experiments.
Returns the final state of the beads and associated bond directions and lattice connections.
If show_E is true, additionally returns the energy over time.


"""
function burnin(conf::Union{RotationConfig, PatchConfig}, Nt::Int; show_E=false)
    conf_burnin = deepcopy(conf)
    conf_burnin = @set conf_burnin.external_force = NoExternalForce()
    beads, bead_info = initialise(conf_burnin)

    step = 20
    time = 0:step:Nt
    E = zeros((6,Nt÷step+1))

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf_burnin, conf_burnin.iter_pars)
        # if i % step == 0 && show_E==true
        #     E[:,i÷step+1] = total_energy(beads, bead_info, dirs, conf_burnin.spring_consts)
        # end
    end

    # if show_E==true
    #     beads, bead_info, dirs, E, time
    # end

    return beads, bead_info
end
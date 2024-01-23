"""
    burnin(conf::Union{MicrotubuleConfig, PatchConfig}, Nt::Int; show_E=false)

Initialise and run simulation for configuration `conf` with `Nt` timesteps. 
Finds equilibrium position to then be used for experiments.
Returns the final state of the beads and associated bond directions and lattice connections.
If show_E is true, additionally returns the energy over time.


"""
function burnin(conf::Union{MicrotubuleConfig, PatchConfig}, Nt::Int)
    conf_burnin = deepcopy(conf)
    conf_burnin = @set conf_burnin.external_force = NoExternalForce()
    beads, bead_info = initialise(conf_burnin)

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf_burnin, conf_burnin.iter_pars)
    end

    return beads, bead_info
end
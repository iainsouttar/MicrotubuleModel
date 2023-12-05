function iterate!(beads, bead_info, conf, dirs)
    update!(beads, bead_info, dirs, conf, conf.iter_pars)
end

function update!(beads, bead_info, dirs, conf, iter::RK4Pars)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = lastindex(beads)
    F = zeros(Float64, (3, lastindex(beads)))
    torque = similar(F)
    b = deepcopy(beads)
    
    k1x, k1q = forward(b, bead_info, F, torque, dirs, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k1x[i]*dt/2
        b[i].q = beads[i].q + (i>N)*k1q[i]*dt/2
    end
    k2x, k2q = forward(b, bead_info, F, torque, dirs, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k2x[i]*dt/2
        b[i].q = beads[i].q + (i>N)*k2q[i]*dt/2
    end
    k3x, k3q = forward(b, bead_info, F, torque, dirs, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k3x[i]*dt
        b[i].q = beads[i].q + (i>N)*k3q[i]*dt
    end
    k4x, k4q = forward(b, bead_info, F, torque, dirs, conf)

    k_x = @. (k1x + 2*k2x + 2*k3x + k4x)*dt/6
    k_q = @. (k1q + 2*k2q + 2*k3q + k4q)*dt/6
    @threads for i in 1:Ntot
        beads[i].x = beads[i].x + (i>N)*k_x[i]
        beads[i].q = beads[i].q + (i>N)*k_q[i]
    end
end

function update!(beads, bead_info, dirs, conf, iter::EulerPars)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    F = zeros(Float64, (3, lastindex(beads)))
    torque = similar(F)
    
    Δx, Δq = forward(beads, bead_info, F, torque, dirs, conf)

    @threads for i in 1:lastindex(beads)
        beads[i].x += (i>N)*Δx[i]*dt
        beads[i].q += (i>N)*Δq[i]*dt
    end
end


function forward(beads, bead_info, F, torque, dirs, conf)
    Ntot = lastindex(beads)
    x = Vector{BeadPos}(undef, Ntot)
    q = Vector{Quaternions.Quaternion}(undef, Ntot)
    @unpack damp_x, damp_theta, dt = conf.iter_pars

    eval_forces_and_torques!(F, torque, beads, bead_info, dirs, conf.spring_consts)

    @tturbo for i in 1:Ntot
        x[i], q[i] = forward_(beads[i], F[:,i], torque[:,i], damp_x, damp_theta)
    end
    return x, q
end


@inline function eval_forces_and_torques!(F, torque, beads, bead_info, dirs, consts)
    K = consts.K
    @threads for i in 1:lastindex(beads)
        F[:,i] = linear_forces(beads[i], beads, bead_info[i], consts)
        torque[:,i] = angular_forces!(F, beads[i], beads, bead_info[i], dirs[bead_info[i].α], K)
    end
end

@inline function forward_(b, F, τ, γ_x, γ_θ)
    #x = - F * γ_x
    
    q = sign(b.q)
    q_τ = quat(0, τ...)
    spin = 0.5 * q * q_τ 
    #q = q + spin * γ_θ
    return - F * γ_x, spin * γ_θ
end


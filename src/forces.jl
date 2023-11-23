

function linear_forces(b1, lattice, consts::SpringConst)
    @unpack k_lat, k_long, k_in, k_in_kin = consts
    @unpack l0_lat, l0_long, l0_in, l0_in_kin = consts
    lat, long = b1.lat_nn, b1.long_nn
    intra = b1.intra_nn
    F = MVector{3,Float64}(0,0,0)
    F += sum(spring_force(b1.x - lattice[b].x, l0_lat, k_lat) for b in lat if b != 0)

    if long != 0
        F += spring_force(b1.x - lattice[long].x, l0_long, k_long)
    end
    if intra != 0
        k, l0 = lattice[intra].kinesin == true ? (k_in_kin, l0_in_kin) : (k_in, l0_in)
        F += spring_force(b1.x - lattice[intra].x, l0, k)
    end

    return F
end

function spring_force(r::BeadPos, l0::Real, k::Real)
    d = norm(r)
    @fastmath ΔL = d - l0
    return (k*ΔL/d) .* r
end

function eval_forces_and_toques!(F, torque, lattice, dirs, consts)
    M = lastindex(lattice)
    K = consts.K
    s = [0.0,0.0,0.0]
    for i in 1:M
        F[:,i] = linear_forces(lattice[i], lattice, consts)
        s .+= F[:,i]
        torque[:,i] = angular_forces!(F, lattice[i], lattice, dirs[lattice[i].α], K)
    end
end

function update_pos_orientation!(lattice, F, torque, pars, N)
    @unpack damp_x, damp_theta, dt, mass, moment_inertia = pars
    dt_x = dt / (damp_x * mass)
    dt_θ = dt / (2 * damp_theta * moment_inertia)
    @threads for i in 1:lastindex(lattice)
        if i > N
            lattice[i].x .-= F[:,i] .* dt_x
            q_τ = quat(0, torque[:,i]...)
            lattice[i].q = sign(lattice[i].q)
            lattice[i].q = lattice[i].q + q_τ * lattice[i].q .* dt_θ
        end
        lattice[i].q = sign(lattice[i].q)
    end
end

function iterate!(lattice, conf, dirs)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N

    F = zeros(Float64, (3, lastindex(lattice)))
    torque = similar(F)

    eval_forces_and_toques!(F, torque, lattice, dirs, conf.spring_consts)

    update_pos_orientation!(lattice, F, torque, conf.iter_pars, N)
end

# going from the bottom left, spiral upwards taking left and below forces
# calculate up and right forces and store 
# each bond has three forces associated
# should know if mutual force has been calculated based on lattice structure
# do a global pass to sum up forces 



spring_energy(r, l0, k) = k/2*(norm(r)-l0)^2

function bead_energy(b1, lattice, consts)
    @unpack k_lat, k_long, k_in, l0_lat, l0_long, l0_in = consts
    lat = b1.lat_nn
    long = b1.long_nn
    intra = b1.intra_nn
    E = 0.0
    
    E += sum(spring_energy(b1.x - lattice[b].x, l0_long, k_long) for b in long if b != 0)
    E += sum(spring_energy(b1.x - lattice[b].x, l0_lat, k_lat) for b in lat  if b != 0)
    E += long == 0 ? 0.0 : spring_energy(b1.x - lattice[long].x, l0_long, k_long)
    E += intra == 0 ? 0.0 : spring_energy(b1.x - lattice[intra].x, l0_in, k_in)

    return E/2
end

function total_energy(lattice, consts::SpringConst)
    return sum(bead_energy(b, lattice, consts) for (idx,b) in enumerate(lattice))
end
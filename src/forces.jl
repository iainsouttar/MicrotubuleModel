
"""
    linear_forces(b1::Bead, lattice::Vector{Bead}, consts::SpringConst)::MVector{3,Float64}

Calculate 3D force acting on bead `b1` due to the linear spring attaching its neighbours.
"""
function linear_forces(b1::Bead, lattice::Vector{Bead}, consts::SpringConst)::MVector{3,Float64}
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

"""
    spring_force(r::BeadPos, l0::Real, k::Real)

Return force due to displacement r of a spring

# Arguments
- `r::BeadPos`: 3D displacement vector
- `l0::Real`: rest length of spring
- `k::Real`: spring stiffness

# Returns
- `MVector{3, Float64}`: directed force
"""
@inline function spring_force(r::BeadPos, l0::Real, k::Real)
    @fastmath d = norm(r)
    @fastmath ΔL = d - l0
    return @fastmath (k*ΔL/d) .* r
end


@inline function eval_forces_and_torques!(F, torque, lattice, dirs, consts)
    K = consts.K
    @threads for i in 1:lastindex(lattice)
        F[:,i] = linear_forces(lattice[i], lattice, consts)
        torque[:,i] = angular_forces!(F, lattice[i], lattice, dirs[lattice[i].α], K)
    end
end


@inline function update_pos_orientation!(lattice, F, torque, pars, N)
    @unpack damp_x, damp_theta, dt, mass, moment_inertia = pars
    dt_x = dt / (damp_x * mass)
    dt_θ = dt / (damp_theta * moment_inertia)
    @threads for i in 1:lastindex(lattice)
        if i > N
            @fastmath step!(lattice[i], F[:,i], torque[:,i], dt_x, dt_θ)
        end
    end
end

@inline function step!(b, F, τ, dt_x, dt_θ)
    b.x .-= F * dt_x
    
    q = sign(b.q)
    q_τ = quat(0, τ...)
    #spin = 0.5* q_τ * q
    spin = 0.5 * q * q_τ
    b.q = q + spin * dt_θ
    return spin
end


function iterate!(lattice, conf, dirs)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N

    F = zeros(Float64, (3, lastindex(lattice)))
    torque = similar(F)

    eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

    update_pos_orientation!(lattice, F, torque, conf.iter_pars, N)
end

# going from the bottom left, spiral upwards taking left and below forces
# calculate up and right forces and store 
# each bond has three forces associated
# should know if mutual force has been calculated based on lattice structure
# do a global pass to sum up forces 



@inline spring_energy(r, l0, k) = k/2*(norm(r)-l0)^2

function bead_energy(b1, lattice, dirs, consts)
    @unpack k_lat, k_long, k_in, l0_lat, l0_long, l0_in = consts
    @unpack l0_in_kin, k_in_kin = consts
    lat = b1.lat_nn
    long = b1.long_nn
    intra = b1.intra_nn
    (k_intra, l0_intra) = b1.kinesin ? (k_in_kin, l0_in_kin) : (k_in, k_in_kin)
    E = 0.0
    
    E += sum(spring_energy(b1.x - lattice[b].x, l0_lat, k_lat) for b in lat  if b != 0)
    E += long == 0 ? 0.0 : spring_energy(b1.x - lattice[long].x, l0_long, k_long)
    E += intra == 0 ? 0.0 : spring_energy(b1.x - lattice[intra].x, l0_intra, k_intra)

    bonds = [intra,lat[2],long, lat[1]]
    for (bond, dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            r = lattice[bond].x - b1.x
            v = orientate_vector(dir, b1.q)
            theta, _, _ = bond_angle(v, r)
            E += consts.K*theta^2
        end
    end

    return E/2
end

function total_energy(lattice, dirs, consts::SpringConst)
    return sum(bead_energy(b, lattice, dirs[b.α], consts) for (idx,b) in enumerate(lattice))
end
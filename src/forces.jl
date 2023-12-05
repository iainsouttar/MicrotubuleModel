
"""
    linear_forces(b1::Bead, lattice::Vector{Bead}, consts::SpringConst)::MVector{3,Float64}

Calculate 3D force acting on bead `b1` due to the linear spring attaching its neighbours.
"""
function linear_forces(b1::Bead, lattice::Vector{Bead}, neighbours::BeadPars, consts::SpringConst)::MVector{3,Float64}
    @unpack k_lat, k_long, k_in, k_in_kin = consts
    @unpack l0_lat, l0_long, l0_in, l0_in_kin = consts
    @unpack north, east, south, west = neighbours
    F = MVector{3,Float64}(0,0,0)
    if north != 0
        F += spring_force(b1.x - lattice[north].x, l0_long, k_long)
    end
    if east != 0
        F += spring_force(b1.x - lattice[east].x, l0_lat, k_lat)
    end
    if south != 0
        k, l0 = b1.kinesin == true ? (k_in_kin, l0_in_kin) : (k_in, l0_in)
        F += spring_force(b1.x - lattice[south].x, l0, k)
    end
    if west != 0
        F += spring_force(b1.x - lattice[west].x, l0_lat, k_lat)
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


# going from the bottom left, spiral upwards taking left and below forces
# calculate up and right forces and store 
# each bond has three forces associated
# should know if mutual force has been calculated based on lattice structure
# do a global pass to sum up forces 



@inline spring_energy(r, l0, k) = k/2*(norm(r)-l0)^2

function bead_energy(b1, beads, neighbours, dirs, consts)
    @unpack k_lat, k_long, k_in, l0_lat, l0_long, l0_in = consts
    @unpack l0_in_kin, k_in_kin = consts
    @unpack north, east, south, west = neighbours
    (k_intra, l0_intra) = b1.kinesin ? (k_in_kin, l0_in_kin) : (k_in, k_in_kin)
    E = 0.0
    
    #E += sum(spring_energy(b1.x - bead[b].x, l0_lat, k_lat) for b in lat  if b != 0)
    E += north == 0 ? 0.0 : spring_energy(b1.x - beads[north].x, l0_long, k_long)
    E += east == 0 ? 0.0 : spring_energy(b1.x - beads[east].x, l0_lat, k_lat)
    E += south == 0 ? 0.0 : spring_energy(b1.x - beads[south].x, l0_intra, k_intra)
    E += west == 0 ? 0.0 : spring_energy(b1.x - beads[west].x, l0_lat, k_lat)

    bonds = [north, east, south, west]
    for (bond, dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            r = beads[bond].x - b1.x
            v = orientate_vector(dir, b1.q)
            theta, _, _ = bond_angle(v, r)
            E += consts.K*theta^2
        end
    end

    return E/2
end

function total_energy(beads, bead_info, dirs, consts::SpringConst)
    return sum(bead_energy(b, beads, b_, dirs[b_.α], consts) for (b,b_) in zip(beads, bead_info))
end
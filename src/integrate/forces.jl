
"""
    linear_forces(b1::Bead, lattice::Vector{Bead}, b_info::BeadPars, consts::SpringConst)::MVector{3,Float64}

Calculate 3D force acting on bead `b1` due to the linear spring attaching its neighbours.
"""
function linear_forces(
    b1::Bead,
    lattice::Vector{Bead},
    b_info::BeadPars, 
    consts::SpringConst
)
    @unpack k_lat, k_long, k_in, k_in_kin = consts
    @unpack l0_lat, l0_long, l0_in, l0_in_kin = consts
    @unpack north, east, south, west = b_info
    F = MVector{3,Float64}(0,0,0)

    k, l0 = b1.kinesin == true ? (k_in_kin, l0_in_kin) : (k_in, l0_in)
    (k_north, l0_north) = b_info.α ? (k_long, l0_long) : (k, l0)
    (k_south, l0_south) = b_info.α ? (k, l0) : (k_long, l0_long)

    if north != 0
        F += spring_force(lattice[north].x - b1.x, l0_north, k_north)
    end
    if east != 0
        F += spring_force(lattice[east].x - b1.x, l0_lat, k_lat)
    end
    if south != 0
        F += spring_force(lattice[south].x - b1.x, l0_south, k_south)
    end
    if west != 0
        F += spring_force(lattice[west].x - b1.x, l0_lat, k_lat)
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
@inline @fastmath function spring_force(r::BeadPos, l0::Real, k::Real)
    d = norm(r)
    ΔL = d - l0
    return (k*ΔL/d) * r
end


# going from the bottom left, spiral upwards taking left and below forces
# calculate up and right forces and store 
# each bond has three forces associated
# should know if mutual force has been calculated based on lattice structure
# do a global pass to sum up forces 

"""
    linear_forces(b1::Bead, lattice::Vector{Bead}, b_info::BeadPars, consts::SpringConst)::MVector{3,Float64}

Calculate 3D force acting on bead `b1` due to the linear spring attaching its neighbours.
"""
function linear_forces(
    b1::Bead,
    bonds,
    b_info::BeadPars
)
    F = MVector{3,Float64}(0,0,0)
    # loop over the bond  parameters and neighbour bead positions
    for (k, l0, b) in zip(b_info.consts, b_info.lengths, bonds)
        F += spring_force(b.x - b1.x, l0, k)
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

@option "spring consts" struct SpringConst
    K_alpha_long::Float64 = 8.5
    K_beta_long::Float64 = 8.5
    K_lat::Float64 = 8.5
    K_alpha_in::Float64 = 8.5
    K_beta_in::Float64 = 8.5
    K_in_kin::Float64 = 0.0
    k_long::Float64 = 3.5
    k_lat::Float64 = 14.0
    k_in::Float64 = 3.5
    k_in_kin::Float64 = 8.0
    l0_long::Float64 = 4.05
    l0_lat::Float64 = 5.21
    l0_in::Float64 = 4.05
    l0_in_kin::Float64 = 4.5
    K_torque_long::Float64 = 1.0
    K_torque_lat::Float64 = 1.0
    K_torque_in::Float64 = 1.0
    K_alpha::Float64 = 10.0
    K_beta::Float64 = 0.0
end

"""
    sort_spring_consts(consts::SpringConst, α::Bool)

Create vectors for spring constants and natural lengths organised by directions: north, east, south, west.
Uses α flag to determine whether north or south bound is intra-dimer bond
"""
function sort_spring_consts(consts::SpringConst, α::Bool)
    @unpack k_lat, k_long, k_in = consts
    @unpack l0_lat, l0_long, l0_in = consts
    @unpack K_alpha_long, K_beta_long, K_lat, K_alpha_in, K_beta_in = consts
    @unpack K_torque_long, K_torque_lat, K_torque_in = consts
    @unpack K_alpha, K_beta = consts




    # flip intrinsic/longitudinal bond based on alpha/beta
    (k_north, l0_north) = α ? (k_long, l0_long) : (k_in, l0_in)
    (k_south, l0_south) = α ? (k_in, l0_in) : (k_long, l0_long)
    (K_north, K_south) = α ? (K_alpha_long, K_alpha_in) : (K_beta_in, K_beta_long)
    (K_torque_north, K_torque_south) = α ? (K_torque_long, K_torque_in) : (K_torque_in, K_torque_long)


    k = [k_north, k_lat, k_south, k_lat]
    l0 = [l0_north, l0_lat, l0_south, l0_lat]
    K = [K_north, K_lat, K_south, K_lat]
    K_torque = [K_torque_north, K_torque_lat, K_torque_south, K_torque_lat]
    K_pos = [K_alpha, K_beta]

    return k, l0, K, K_torque, K_pos
end

###################################################################

struct BeadPars
    α::Bool
    bonds::Vector{Int}
    directions::Vector{SVector{3,Float64}}
    lin_consts::Vector{Float64}
    bend_consts::Vector{Float64}
    lengths::Vector{Float64}
    torque_consts::Vector{Float64}
    pos_consts::Float64
end

mutable struct Lattice{N}
    x::MVector{N, BeadPos}
    q::Vector{Quaternions.Quaternion}
    kinesin::MMatrix{N, 2, Bool}
    rates::MVector{N, Float64}
end

function Lattice(x, q, kinesin, rates)
    N = length(rates)
    print(N)
    return Lattice{N}(x, q, kinesin, rates)
end


Base.size(l::Lattice) = size(l.kinesin)

Base.length(l::Lattice) = length(l.kinesin)


"""
    create_patch(dirs, consts, N_lat::Int, N_long::Int, a::Real, δx::Real; S::Int=3, N::Int=13)

Construct a lattice of beads connected by springs.

# Arguments
- `dirs::Dict`: bond directions for either alpha/beta tubulin
- `consts::SpringConst`: spring constants
- `N_lat::Int`: number of beads laterally
- `N_long::Int`: number of beads in protofilaments
- `a::Real`: intial separation along protofilament
- `δx::Real`: initial separation laterally

# Optional Arguments
- `S::Int=3`: helix-start number (longitudinal increase between first and last in a ring)
- `N::Int=13`: protofilament number (number of beads in each ring)

# Returns
- `Lattice`: lattice of connected beads
- `Vector{BeadPars}`: bead connections and alpha/beta type
"""
function create_patch(
    dirs::Dict,
    consts::SpringConst,
    N_lat::Int,
    N_long::Int,
    a::Real,
    δx::Real;
    S::Int=3,
    N::Int=13
)
    bead_info = Vector{BeadPars}(undef,N_lat*N_long)

    r = S*a/N  # subunit rise
    R = N*δx/2π  # radius of cylinder

    angles = repeat(range(0.0,2π*(N-1)/N,N)[1:N_lat],N_long)
    vertical_offset = repeat(0:N_long-1,inner=N_lat)

    x = [set_pos(θ, idx, R, r, z_0=a*z, N=N_lat) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    q = Vector{Quaternions.Quaternion}(undef, N_lat*N_long)
    kinesin = MMatrix{N_lat*N_long,2,Bool}(zeros(2*N_long*N_lat))
    alpha = reshape([Bool(z%2==1) for z in vertical_offset], (N_lat,N_long))
    rates = MVector{N_lat*N_long, Float64}(zeros(N_lat*N_long))
    #the 0 in the reshape used to be a 1
    print(size(kinesin))

    for i in 1:N_long
        for j in 1:N_lat
            idx = (i-1)*N_lat+j
            if N_lat==N
                lat, (north, south) = neighbours(idx, N_long*N_lat, N=N, S=S)
            else
                lat = lateral_nn_patch(idx, N_lat)
                (north, south) = long_nn_patch(idx, N_lat, N_long)
            end

            q[idx] = quat_from_axisangle([0,0,1],-π/2+angles[idx]) #+(i-1)*π/4)
            # bonds organised by north, east, south, west
            bonds = [north, lat[2], south, lat[1]]

            k, l0, K, K_torque, K_pos = sort_spring_consts(consts, alpha[j,i])
            # only attach the bead to beads which exist
            # number of bonds not always 4!
            nonzero = bonds .!= 0


            K_position = alpha[j,i] ? K_pos[1] : K_pos[2]
            # construct all info about bond stiffness & directions etc
            bead_info[idx] = BeadPars(
                alpha[j,i], 
                unique(bonds[nonzero]),
                dirs[alpha[j,i]][nonzero],
                k[nonzero],
                K[nonzero],
                l0[nonzero],
                K_torque[nonzero],
                K_position
            )
        end

        #q[1] = quat_from_axisangle([0,0,1], )
    end
    #l = Lattice(x, q, kinesin, rates)

    return Lattice(x, q, kinesin, rates), bead_info
end


"""
    create_lattice(dirs, consts, num_rings::Int, a::Real, δx::Real; S=3, N=13)

Construct a full lattice of beads connected by springs.

# Arguments
- `dirs::Dict`: bond directions for either alpha/beta tubulin
- `consts::SpringConst`: spring constants
- `num_rings::Int`: number of beads laterally
- `a::Real`: intial separation along protofilament
- `δx::Real`: initial separation laterally

# Returns
- `Lattice`: lattice of connected beads
- `Vector{BeadPars}`: bead connections and alpha/beta type
"""
create_lattice(dirs, consts, num_rings, a, δx; S::Int=3, N::Int=13) = create_patch(dirs, consts, N, num_rings, a, δx; S=S, N=N)


"""
    create_dimer(dirs, consts, a::Real)

Construct a full lattice of beads connected by springs.

# Arguments
- `dirs::Dict`: bond directions for either alpha/beta tubulin
- `consts::SpringConst`: spring constants
- `a::Real`: intial separation along protofilament

# Returns
- `Lattice`: lattice of the two connected beads
- `Vector{BeadPars}`: bead connections and alpha/beta type
"""
function create_dimer(dirs, consts, a::Real)
    x = [BeadPos(0,0,0), BeadPos(0,0,a)]
    q = quat_from_axisangle([0,0,1],-π/2)
    kinesin = [false, false]
    @unpack k_in, l0_in, K_in = consts

    bead_info = [
        BeadPars(false, [2],[dirs[false][1]],[k_in],[K_in], [l0_in]),
        BeadPars(true, [1], [dirs[true][3]], [k_in], [K_in], [l0_in])
    ]

    return Lattice(x, [q,q], kinesin), bead_info
end


"""
    create_PF(num_rings::Int, a::Real, δx::Real; S=3, N=13)

Construct a single protofilament of beads connected by springs.

# Arguments
- `dirs::Dict`: bond directions for either alpha/beta tubulin
- `consts::SpringConst`: spring constants
- `num_rings::Int`: number of beads laterally
- `a::Real`: intial separation along protofilament
- `δx::Real`: initial separation laterally

# Returns
- `Lattice`: lattice of connected beads
- `Vector{BeadPars}`: bead connections and alpha/beta type
"""
function create_PF(dirs, consts, n::Int, a::Real)
    bead_info = Vector{BeadPars}(undef,n)

    x = [BeadPos(0,0,i*a) for i in 0:n-1]
    q = quat_from_axisangle([0,0,1],-π/2)
    kinesin = [false for i in 1:n]
    @unpack k_in, k_long, l0_in, l0_long, K_in, K_long = consts

    bead_info[1] = BeadPars(false, [2],[dirs[false][1]],[k_in],[K_in], [l0_in])
    for i in 2:n-1
        α = Bool(i%2==0)
        bead_info[i] = BeadPars(
            α, 
            [i+1,i-1],
            [dirs[α][1],dirs[α][3]],
            α ? [k_long, k_in] : [k_in, k_long],
            α ? [K_long, K_in] : [K_in, K_long],
            α ? [l0_long, l0_in] : [l0_in, l0_long],
        )
    end
    bead_info[n] = BeadPars(true, [n-1],[dirs[true][3]],[k_in],[K_in], [l0_in])

    return Lattice(x, [q for i in 1:n], kinesin), bead_info
end


"""
    set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos

Set position of the bead in a lattice.
"""
function set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos
    return BeadPos(
        R*cos(θ_0-θ),
        R*sin(θ_0-θ), #+(idx)/4,
        z_0 + mod(idx-1,N)*r
    )
end

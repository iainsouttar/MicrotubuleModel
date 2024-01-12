"""
    set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos

Set position of the bead in a lattice.
"""
function set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos
    return BeadPos(
        R*cos(θ_0-θ),
        R*sin(θ_0-θ),
        z_0 + mod(idx-1,N)*r
    )
end

"""
    create_patch(N_lat::Int, N_long::Int, a::Real, δx::Real; S::Int=3, N::Int=13)

Construct a lattice of beads connected by springs.

# Arguments
- `N_lat::Int`: number of beads laterally
- `N_long::Int`: number of beads in protofilaments
- `a::Real`: intial separation along protofilament
- `δx::Real`: initial separation laterally

# Returns
- `Vector{Bead}`: lattice of connected beads
"""
function create_patch(dirs, consts, N_lat::Int, N_long::Int, a::Real, δx::Real; S::Int=3, N::Int=13)
    #beads = Vector{Bead}(undef,N_lat*N_long)
    #lattice = Lattice{Int(N_lat*N_long)}()
    bead_info = Vector{BeadPars}(undef,N_lat*N_long)

    r = S*a/N
    R = N*δx/2π

    angles = repeat(range(0.0,2π*(N-1)/N,N)[1:N_lat],N_long)
    vertical_offset = repeat(0:N_long-1,inner=N_lat)

    x = [set_pos(θ, idx, R, r, z_0=a*z, N=N_lat) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    q = Vector{Quaternions.Quaternion}(undef, N_lat*N_long)
    kinesin = MVector{N_lat*N_long, Bool}(zeros(Bool,N_lat*N_long))
    alpha = reshape([Bool(z%2==1) for z in vertical_offset], (N_lat,N_long))

    for i in 1:N_long
        for j in 1:N_lat
            idx = (i-1)*N_lat+j
            if N_lat==N
                lat, (north, south) = neighbours(idx, N_long*N_lat)
            else
                lat = lateral_nn_patch(idx, N_lat)
                (north, south) = long_nn_patch(idx, N_lat, N_long)
            end

            q[idx] = quat_from_axisangle([0,0,1],-π/2+angles[idx])
            #beads[idx] = Bead(x[idx], q, false)
            #q[idx] = q

            bonds = [north, lat[2], south, lat[1]]

            @unpack k_lat, k_long, k_in, k_in_kin = consts
            @unpack l0_lat, l0_long, l0_in, l0_in_kin = consts

            (k_north, l0_north) = alpha[j,i] ? (k_long, l0_long) : (k_in, l0_in)
            (k_south, l0_south) = alpha[j,i] ? (k_in, l0_in) : (k_long, l0_long)

            k = [k_north, k_lat, k_south, k_lat]
            l0 = [l0_north, l0_lat, l0_south, l0_lat]

            bead_info[idx] = BeadPars(
                alpha[j,i], 
                bonds[bonds .!= 0],
                dirs[alpha[j,i]][bonds .!= 0],
                k[bonds .!= 0],
                l0[bonds .!= 0]
            )
        end
    end


    return Lattice(x, q, kinesin), bead_info
end


"""
    create_lattice(num_rings::Int, a::Real, δx::Real; S=3, N=13)

Construct a full lattice of beads connected by springs.

# Arguments
- `num_rings::Int`: number of beads laterally
- `a::Real`: intial separation along protofilament
- `δx::Real`: initial separation laterally

# Returns
- `Vector{Bead}`: lattice of connected beads
"""
create_lattice(dirs, consts, num_rings, a, δx; S::Int=3, N::Int=13) = create_patch(dirs, consts, N, num_rings, a, δx; S=S, N=N)


function create_dimer(dirs, consts, a::Real, δx::Real; S::Int=3, N::Int=13)
    bead_info = Vector{BeadPars}(undef,2)

    r = S*a/N
    R = N*δx/2π

    angles = [0.0,0.0]#range(0.0,2π*(N-1)/N,N)[1:13:14]
    vertical_offset = [0,1]

    x = [BeadPos(0,0,0), BeadPos(0,0,a)]
    #x = [set_pos(θ, idx, R, r, z_0=a*z, N=N) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    q = quat_from_axisangle([0,0,1],-π/2)
    kinesin = [false, false]

    @unpack k_in, l0_in = consts

    bead_info[1] = BeadPars(
        false, 
        [2],
        [dirs[false][1]],
        [k_in],
        [l0_in]
    )
    bead_info[2] = BeadPars(
        true, 
        [1],
        [dirs[true][3]],
        [k_in],
        [l0_in]
    )
    return Lattice(x, [q,q], kinesin), bead_info
end



function create_PF(dirs, consts, n::Int, a::Real, δx::Real; S::Int=3, N::Int=13)
    bead_info = Vector{BeadPars}(undef,n)

    x = [BeadPos(0,0,i*a) for i in 0:n-1]

    q = quat_from_axisangle([0,0,1],-π/2)
    kinesin = [false for i in 1:n]

    @unpack k_in, l0_in = consts

    bead_info[1] = BeadPars(
        false, 
        [2],
        [dirs[false][1]],
        [k_in],
        [l0_in]
    )
    for i in 2:n-1
        bead_info[i] = BeadPars(
            Bool(i%2==0), 
            [i+1,i-1],
            [dirs[Bool(i%2==0)][1],dirs[Bool(i%2==0)][3]],
            [k_in, k_in],
            [l0_in, l0_in]
        )
    end
    bead_info[n] = BeadPars(
        true, 
        [n-1],
        [dirs[true][3]],
        [k_in],
        [l0_in]
    )
    return Lattice(x, [q for i in 1:n], kinesin), bead_info
end
"""
    set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos

Set position of the bead in a lattice.
"""
function set_pos(θ::Real, idx::Int, R::Real, r::Real; θ_0=0.0, z_0=0.0, N=13)::BeadPos
    return BeadPos(
        R*cos(θ_0+θ),
        R*sin(θ_0+θ),
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
function create_patch(N_lat::Int, N_long::Int, a::Real, δx::Real; S::Int=3, N::Int=13)
    beads = Vector{Bead}(undef,N_lat*N_long)
    bead_info = Vector{BeadPars}(undef,N_lat*N_long)

    r = S*a/N
    R = N*δx/2π

    angles = repeat(range(0.0,2π*(N-1)/N,N)[1:N_lat],N_long)
    vertical_offset = repeat(0:N_long-1,inner=N_lat)

    x = [set_pos(θ, idx, R, r, z_0=a*z, N=N_lat) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
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

            q = quat_from_axisangle([0,0,1],-π/2-angles[idx])
            beads[idx] = Bead(x[idx], q, false)

            bead_info[idx] = BeadPars(
                alpha[j,i], 
                north, lat[2], south, lat[1]
            )
        end
    end
    return beads, bead_info
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
create_lattice(num_rings, a, δx; S::Int=3, N::Int=13) = create_patch(N, num_rings, a, δx; S=S, N=N)


function create_dimer(a::Real, δx::Real; S::Int=3, N::Int=13)
    beads = Vector{Bead}(undef,2)
    bead_info = Vector{BeadPars}(undef,2)

    r = S*a/N
    R = N*δx/2π

    angles = range(0.0,2π*(N-1)/N,N)[1:2]
    vertical_offset = [0,0]

    x = [set_pos(θ, idx, R, r, z_0=a*z, N=N) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    q = quat_from_axisangle([0,0,1],-π/2)
    # q = quat_from_axisangle([0,0,1],-π/2-angles[idx])

    beads[1] = Bead(x[1], copy(q), false)
    beads[2] = Bead(x[2], copy(q), false)

    bead_info[1] = BeadPars(
        false, 
        0, 2, 0, 0
    )
    bead_info[2] = BeadPars(
        false, 
        0, 0, 0, 1
    )
    return beads, bead_info
end


"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the lateral direction for bead `idx` of `total`
"""
function lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}
    if idx % 13 == 0 
        if idx > total-2*13-1
            return (idx-1, 0)
        end
        return (idx-1, idx+1+2*13)
    elseif idx % 13 == 1
        if idx < 2*13+1
            return (0, idx+1)
        end
        return (idx-1-2*13, idx+1)
    end
    return (idx-1, idx+1)
end

"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the longitudinal direction for bead `idx` of `total`
"""
function long_nn(idx::Int, total::Int)::Tuple{Int,Int}
    if idx<14
        return idx+13, 0
    elseif idx>total-13
        return 0, idx-13
    end
    return idx+13, idx-13
end

"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the lateral direction for bead `idx` of `total` in a patch
"""
function lateral_nn_patch(idx::Int, N_lat::Int)::Tuple{Int,Int}
    if idx % N_lat == 0 
        return (idx-1, 0)
    elseif idx % N_lat == 1
        return (0, idx+1)
    end
    return (idx-1, idx+1)
end


"""
    long_nn_patch(idx::Int, N_lat::Int, N_long::Int)::Tuple{Int,Int}

Return neighbours in the long direction for bead `idx` of `total` in a patch
"""
function long_nn_patch(idx::Int, N_lat::Int, N_long::Int)::Tuple{Int,Int}
    if idx<N_lat+1
        return idx+N_lat, 0
    elseif idx>N_lat*(N_long-1)
        return 0, idx-N_lat
    end
    return idx+N_lat, idx-N_lat
end


"""
    neighbours(idx::Int, total::Int)

Return neighbours in the lateral and longitudinal direction for bead `idx` of `total`
"""
function neighbours(idx::Int, total::Int)::Tuple{Tuple{Int,Int}, Tuple{Int, Int}}
    return lateral_nn(idx, total), long_nn(idx, total)
end


function set_bond_angles(conf)
    @unpack S, N, dx, a = conf.lattice
    r = S*a/N
    R = N*dx/2π
    l = 2*R*sin(π/N)
    ϕ = atan(r,l)
    #δ = 0.2
    δ = 0.0

    conf = @set conf.alpha = MicrotubuleSpringModel.AlphaConfirm(
        [π/2, -δ],
        [π/13, π/2-ϕ],
        [π/2, π],
        [π-π/13,π/2+ϕ]
    )
    conf = @set conf.beta = MicrotubuleSpringModel.BetaConfirm(
        [π/2, 0],
        [π/13, π/2-ϕ],
        [π/2, π+δ],
        [π-π/13,π/2+ϕ]
    )
    return conf
end
function set_pos(θ, idx, R, r; θ_0=0.0, z_0=0.0, N=13)
    return BeadPos(
        R*cos(θ_0+θ),
        R*sin(θ_0+θ),
        z_0 + mod(idx-1,N)*r
    )
end

function create_patch(N_lat, N_long, a, δx; S=3, N=13)
    beads = Vector{Bead}(undef,N_lat*N_long)

    r = S*a/N
    R = N*δx/2π

    angles = repeat(range(0.0,2π*(N-1)/N,N)[1:N_lat],N_long)
    vertical_offset = repeat(0:N_long-1,inner=N_lat)

    positions = [set_pos(θ, idx, R, r, z_0=a*z, N=N_lat) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    alpha = reshape([Bool(z%2==1) for z in vertical_offset], (N_lat,N_long))

    for i in 1:N_long
        for j in 1:N_lat
            if N_lat==N
                lat, (intra, long) = neighbours((i-1)*N_lat+j, N_long*N_lat)
            else
                lat = lateral_nn_patch((i-1)*N_lat+j, N_lat)
                (intra, long) = long_nn_patch((i-1)*N_lat+j, N_lat, N_long)
            end
            long, intra = alpha[j,i]==true ? (long, intra) : (intra, long)

            q = quat_from_axisangle([0,0,1],-π/2-angles[(i-1)*N_lat+j])
            beads[(i-1)*N_lat+j] = Bead(
                positions[(i-1)*N_lat+j], 
                q,
                alpha[j,i], 
                false,
                lat, intra, long
            )
        end
    end
    return beads
end

create_lattice(num_rings, a, δx; S=3, N=13) = create_patch(N, num_rings, a, δx; S=S, N=N)

function lateral_nn(idx, total)::Tuple{Int,Int}
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

function long_nn(idx, total)::Tuple{Int,Int}
    if idx<14
        return idx+13, 0
    elseif idx>total-13
        return 0, idx-13
    end
    return idx+13, idx-13
end


function lateral_nn_patch(idx, N_lat)::Tuple{Int,Int}
    if idx % N_lat == 0 
        return (idx-1, 0)
    elseif idx % N_lat == 1
        return (0, idx+1)
    end
    return (idx-1, idx+1)
end

function long_nn_patch(idx, N_lat, N_long)::Tuple{Int,Int}
    if idx<N_lat+1
        return idx+N_lat, 0
    elseif idx>N_lat*(N_long-1)
        return 0, idx-N_lat
    end
    return idx+N_lat, idx-N_lat
end

function neighbours(idx,total)
    return lateral_nn(idx, total), long_nn(idx, total)
end


function set_bond_angles(conf)
    @unpack S, N, dx, a = conf.lattice
    r = S*a/N
    R = N*dx/2π
    l = 2*R*sin(π/N)
    ϕ = atan(r,l)

    conf = @set conf.alpha = MicrotubuleSpringModel.AlphaConfirm(
        [π/2, -0.2],
        [π/13, π/2-ϕ],
        [π/2, π],
        [π-π/13,π/2+ϕ]
    )
    conf = @set conf.beta = MicrotubuleSpringModel.BetaConfirm(
        [π/2, 0],
        [π/13, π/2-ϕ],
        [π/2, π+0.2],
        [π-π/13,π/2+ϕ]
    )
    return conf
end
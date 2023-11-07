function set_pos(θ, idx, R, r; θ_0=0.0, z_0=0.0)
    return BeadPos(
        R*cos(θ_0+θ),
        R*sin(θ_0+θ),
        z_0 + mod(idx-1,13)*r
    )
end


function create_lattice(num_rings, a, δx; S=3, N=13)
    beads = Vector{Bead}(undef,num_rings*N)

    r = S*a/N
    R = N*δx/2π

    angles = repeat(range(0.0,2π*(N-1)/N,N),num_rings)
    vertical_offset = repeat(0:num_rings-1,inner=N)


    positions = [set_pos(θ, idx, R, r, z_0=z*a) for (idx,(θ,z)) in enumerate(zip(angles,vertical_offset))]
    alpha = reshape([Bool(z%2==0) for z in vertical_offset], (N,num_rings))

    for i in 1:num_rings
        for j in 1:N
            lat, intra, long = neighbours((i-1)*N+j, num_rings*N)
            if alpha==false
                long, intra = intra, long
            end
            beads[(i-1)*N+j] = Bead(
                positions[(i-1)*N+j], 
                BeadAngle(0,0,angles[(i-1)*N+j]), 
                alpha[j,i], 
                lat, intra, long
                )
        end
    end
    return beads
end

function lateral_nn(idx, total)
    if idx % 13 == 0 
        if idx > total-2*13-1
            return (idx-1, 0)
        end
        return (idx-1, idx+1+2*13)
    elseif idx % 13 == 1
        if idx < 2*26+1
            return (0, idx+1)
        end
        return (idx-1-2*13, idx+1)
    end
    
    return (idx-1, idx+1)
end

function long_nn(idx, total)
    if idx<14
        return idx+13, 0
    elseif idx>total-13
        return 0, idx-13
    end
    
    return idx+13, idx-13
end

# return lateral inter bonds, intradimer bond, longitudinal bond
function neighbours(idx,total)
    if idx == 1
        return (idx+1, 0), idx+13, 0
    elseif idx == total
        return (idx-1, 0), 0, idx-13
    end

    if idx<14
        if idx == 13
            return (idx-1, idx+1+2*13), idx+13, 0
        end
        return (idx-1, idx+1), idx+13, 0
    elseif idx>total-13
        if idx == total-12
            return (idx-1-2*13, idx+1), 0, idx-13
        end
        return (idx-1, idx+1), 0, idx-13
    end

    if idx % 13 == 0 
        if idx > total-2*13-1
            return (idx-1, 0), idx+13, idx-13
        end
        return (idx-1, idx+1+2*13), idx+13, idx-13
    elseif idx % 13 == 1
        if idx < 2*26+1
            return (idx+1, 0), idx+13, idx-13
        end
        return (idx-1-2*13, idx+1), idx+13, idx-13
    end
    
    return (idx-1, idx+1), idx+13, idx-13
end
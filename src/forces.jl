
struct LinSpringConst
    k_long::Float64
    k_lat::Float64
    k_intra::Float64
    l0_long::Float64
    l0_lat::Float64
    l0_intra::Float64
end

function linear_forces(b1, lattice, consts::LinSpringConst)
    @unpack k_lat, k_long, k_intra, l0_lat, l0_long, l0_intra = consts
    lat, long = b1.lat_nn, b1.long_nn
    intra = b1.intra_nn
    F = MVector{3,Float64}(0,0,0)
    F += sum(spring_force(b1.x - lattice[b].x, l0_lat, k_lat) for b in lat if b != 0)

    if long != 0
        F += spring_force(b1.x - lattice[long].x, l0_long, k_long)
    end
    if intra != 0
        F += spring_force(b1.x - lattice[intra].x, l0_intra, k_intra)
    end

    return F
end

function spring_force(r::BeadPos, l0::Real, k::Real)
    d = sqrt(dot(r,r))
    @fastmath ΔL = d - l0
    return (k*ΔL/d) .* r
end

function spring_energy(r, l0, k)
    d = sqrt(dot(r,r))
    return k/2*(d-l0)^2
end

function bead_energy(b1, lattice, consts)
    @unpack k_lat, k_long, k_intra, l0_lat, l0_long, l0_intra = consts
    lat = b1.lat_nn
    long = b1.long_nn
    intra = b1.intra_nn
    E = 0.0
    
    E += sum(spring_energy(b1.x - lattice[b].x, l0_long, k_long) for b in long if b != 0)
    E += sum(spring_energy(b1.x - lattice[b].x, l0_lat, k_lat) for b in lat  if b != 0)
    E += long == 0 ? 0.0 : spring_energy(b1.x - lattice[long].x, l0_long, k_long)
    E += intra == 0 ? 0.0 : spring_energy(b1.x - lattice[intra].x, l0_intra, k_intra)

    return E/2
end

function total_energy(lattice, consts::LinSpringConst)
    return sum(bead_energy(b, lattice, consts) for (idx,b) in enumerate(lattice))
end

function iterate!(lattice, consts, dt; S=3, N=13)
    @threads for i in 1:lastindex(lattice)
        if i > N
            F = linear_forces(lattice[i], lattice, consts)
            lattice[i].x .-= F .* dt
        end
    end
end
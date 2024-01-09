
@fastmath bend_energy(r, v, K) = K*(1.0-dot(r,v))

@fastmath spring_energy(r::Real, l0::Real, k::Real) = 0.5*k*(r-l0)^2

"""
    bead_energy(b1::Bead,
        beads::Vector{Bead},
        b_info::BeadPars,
        dirs,
        consts::SpringConst
    )

    Calculate the components of the total energy for the whole lattice.
    Returns E = [lin lateral, lin intrinsic, lin longitudinal, bend lateral, bend intrinsic, bend longitudinal]
"""
function bead_energy(
    b1::Bead,
    beads::Vector{Bead},
    b_info::BeadPars,
    dirs,
    consts::SpringConst
)
    @unpack k_lat, k_long, k_in, k_in_kin = consts
    @unpack l0_lat, l0_long, l0_in, l0_in_kin = consts
    @unpack north, east, south, west = b_info

    k, l0 = b1.kinesin == true ? (k_in_kin, l0_in_kin) : (k_in, l0_in)
    (k_north, l0_north, idx1) = b_info.α ? (k_long, l0_long, 2) : (k, l0, 3)
    (k_south, l0_south, idx2) = b_info.α ? (k, l0, 3) : (k_long, l0_long, 2)

    E = MVector{6,Float64}(zeros(6))

    indices = [idx1, 1, idx2, 1]
    bonds = [north, east, south, west]
    k_ = [k_north, k_lat, k_south, k_lat]
    l0_ = [l0_north, l0_lat, l0_south, l0_lat]

    for (i, bond, dir, k_b, l0_b) in zip(indices, bonds, eachcol(dirs), k_, l0_)
        if bond != 0
            r = beads[bond].x - b1.x
            v = orientate_vector(dir, b1.q)
            E[i] += spring_energy(norm(r), l0_b, k_b)
            E[i+3] += bend_energy(normalize(r), v, consts.K)
        end
    end

    return E ./ 2
end

function total_energy(beads, bead_info, dirs, consts::SpringConst)::Vector{Float64}
    return sum(bead_energy(b, beads, b_, dirs[b_.α], consts) for (b,b_) in zip(beads, bead_info))
end
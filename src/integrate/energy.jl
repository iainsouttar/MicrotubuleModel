
@fastmath bend_energy(r, v, K) = K*(1.0-dot(r,v))

@fastmath spring_energy(r::Real, l0::Real, k::Real) = 0.5*k*(r-l0)^2

"""
    bead_energy(
        x,
        q::Quaternions.Quaternion,
        bonds,
        b::BeadPars,
        K::Real
    )

    Calculate the components of the total energy for the whole lattice.
    Returns E = [lin lateral, lin intrinsic, lin longitudinal, bend lateral, bend intrinsic, bend longitudinal]
"""
function bead_energy(
    x,
    q::Quaternions.Quaternion,
    indices,
    bonds,
    b::BeadPars
)
    E = MVector{6,Float64}(zeros(6))
    #E = 0.0

    for (k, K, l0, dir, bx, idx) in zip(b.lin_consts, b.bend_consts, b.lengths, b.directions, bonds, indices)
        # transform bond direction according to bead orientation
        v = orientate_vector(dir, sign(q))
        @fastmath r = bx - x
        d = norm(r)
        E[idx] += spring_energy(d, l0, k)
        E[idx+3] += bend_energy(r./d, v, K)
    end

    return E ./ 2
end

function total_energy(lattice, bead_info)
    Ntot = length(lattice)
    E = zeros(Float64, (6, Ntot))
    #E = zeros(Float64, length(lattice))
    @inbounds @fastmath @threads for i in 1:length(lattice)
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        indices = bond_indices(i, Ntot, b.α)
        E[:,i] = bead_energy(lattice.x[i], lattice.q[i], indices, bonds, b)
    end
    return vec(sum(E, dims=2))
end
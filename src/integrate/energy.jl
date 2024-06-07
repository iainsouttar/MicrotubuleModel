
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
    b::BeadPars,
    bond_orients
)
    E = MVector{7,Float64}(zeros(7))
    # print(indices, "\n")
    for (k, K, l0, dir, bx, idx, bq, K_torque) in zip(b.lin_consts, b.bend_consts, b.lengths, b.directions, bonds, indices, bond_orients, b.torque_consts)
        # transform bond direction according to bead orientation
        v = orientate_vector(dir, sign(q))
        @fastmath r = bx - x
        d = norm(r)
        E[idx] += spring_energy(d, l0, k)
        E[idx+3] += bend_energy(r./d, v, K)
        rhat, d = norm_and_mag(bx-x)
        e2 = SVector{3,Float64}(1,0,0)
        K_torque = 0
        normal_vec = normalize(cross(rhat, e2))
        orientationi = orientate_vector(e2, sign(q))
        orientationj = orientate_vector(e2, sign(bq))
        orientationiproj = orientationi - dot(orientationi, rhat)*rhat
        orientationjproj = orientationj - dot(orientationj, rhat)*rhat
        #E[7] = norm()
        orientationiproj = orientationiproj/norm(orientationiproj)
        orientationjproj = orientationjproj/norm(orientationjproj)

        #print("\n", orientationiproj,orientationjproj, "\n")
        #rdotv = max(-1, min(1, dot(orientationiproj,orientationjproj))) 
        rdotv = dot(orientationiproj,orientationjproj)

        # τ = K*sin(ϕ)*n̂

        #E[7] = 0.01*K_torque*(acos(clamp(rdotv, -1, 1)))^2
        
        #print("Energy:", E[7]," ", rdotv, " ", x, bx, "\n")
    end
    #print(E[7])


    return E ./ 2
end

function total_energy(lattice, bead_info, N, S)
    Ntot = length(lattice.x)
    E = zeros(Float64, (7, Ntot))
    @inbounds @fastmath @threads for i in 1:length(lattice.x)
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        bond_orients = lattice.q[b.bonds]

        indices = bond_indices(i, Ntot, b.α, N, S)
        E[:,i] = bead_energy(lattice.x[i], lattice.q[i], indices, bonds, b, bond_orients)
    end
    #lngth = orientate_vector((0,0,1), lattice.q[1])
    lngth = norm(lattice.x[1]-lattice.x[Ntot])
    #lngth = norm(lattice.x[1][1])
    return vec(sum(E, dims=2)), 2*E[7, :], lngth[1], E
end
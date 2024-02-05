"""
    internal_forces_and_torques!(F, F_, torque, lattice, bead_info, consts)

Calculate then aggregate the forces and torques on each bead due to internal interactions.
"""
@inline function internal_forces_and_torques!(F, F_, torque, lattice, bead_info)
    N = length(lattice)

    @inbounds @fastmath @threads for i in 1:N
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        bond_orients = lattice.q[b.bonds]
        torque[:,i], F[:,i], F_[:,:,i] = bending_and_spring_forces(lattice.x[i], lattice.q[i], bonds, b, bond_orients)


    end
    @threads for i in 1:N
        b = bead_info[i]
        for j in 1:lastindex(b.bonds)
            @. F[:, b.bonds[j]] += F_[:, j, i]
        end
    end

end



"""
    bending_and_spring_forces(x, q, bonds, b, K)

Calculate 3D torque and force acting on bead `b1` and its neighbours due to the bond angle bending and spring force at `b1`. Updates overall force vectors and returns torque on `b1`.
"""
function bending_and_spring_forces(x, q, bonds, b, bond_orients)
    @unpack lin_consts, bend_consts, lengths, directions = b
    F_ = MMatrix{3, 4, Float64}(undef)
    torque = MVector{3,Float64}(0,0,0)
    F1 = MVector{3,Float64}(0,0,0)
    K_torque = 1
    for (i, (k, K, l0, dir, bx, bq)) in enumerate(zip(lin_consts, bend_consts, lengths, directions, bonds, bond_orients))
        rhat, d = norm_and_mag(bx-x)
        # transform bond direction according to bead orientation
        v = orientate_vector(dir, sign(q))
        # torque from diff between rest direction v and actual r
        τ, F = bending_torque_and_force(v, rhat, d, K)

        #torque from torsion
        orientationi = orientate_vector(BondDirec(1.0,0.0,0.0), sign(q))
        orientationj = orientate_vector(BondDirec(1.0,0.0,0.0), sign(bq))
        τ1 = torsion_torque(rhat, d, orientationi,  orientationj, K_torque)
        
        
        @. F_[:,i] = -F
        torque += τ + τ1
        F1 += F
        F1 += spring_force(rhat, d, l0, k)
    end
    return torque, F1, F_
end


"""
    torque_and_force(v, rhat, d, K)

Calculate the torque and force due to the bond angle bending between the natural direction `v` and the actual direction `rhat` with length `d`.
"""
@inline @fastmath function bending_torque_and_force(v, rhat, d, K)
    if all(rhat .≈ v)
        return zeros(3), zeros(3)
    end

    rdotv = min(1, abs(dot(rhat,v)))
    # τ = K*sin(θ)*n̂
    torque = K*sqrt(1.0-rdotv^2)*normalize(cross(rhat, v))
    # F = -K*sin(θ)θ̂/|r|
    force = -K/d*(v - rdotv*rhat)
    return torque, force
end


"""
    torsion_torque(v, rhat, d, K)

Calculate the torque due to the bond angle torsion between two particles i and j, given the bond direction and the two orientations. Force will always be 0
"""
@inline @fastmath function torsion_torque(rhat, d, zi, zj, K_torque)
    if any(all(dot(zi,zj) .≈ 1))
        return zeros(3)
    end
    

    rdotv = min(1, abs(dot(zi,zj)))
    # τ = K*sin(ϕ)*n̂
    torque = K_torque*sqrt(1.0-rdotv^2)*rhat/d
    return torque
end


function norm_and_mag(r)
    d = norm(r)
    rhat = r ./ d
    return rhat, d
end

"""
    spring_force(rhat, d, l0::Real, k::Real)

Return force due to displacement r of a spring

# Arguments
- `rhat::BeadPos`: 3D direction vector of bond
- `d::Real`: length of bond
- `l0::Real`: rest length of spring
- `k::Real`: spring stiffness

# Returns
- `MVector{3, Float64}`: directed force
"""
@fastmath spring_force(rhat, d::Real, l0::Real, k::Real) = k*(d - l0) * rhat
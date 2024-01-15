function norm_and_mag(r)
    d = norm(r)
    rhat = r ./ d
    return rhat, d
end

"""
    torque_and_force(v, r, K)

Calculate the torque and force due to the bond angle bending between the natural direction v and the actual direction r.
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
    angular_forces!(F, i, b1::Bead, beads, bead_info, dirs, K)::MVector{3,Float64}

Calculate 3D torque and force acting on bead `b1` and its neighbours due to the bond angle bending at `b1`. Updates overall force vectors and returns torque on b1.
"""
function bending_and_spring_forces(x, q, bonds, b, K)
    @unpack consts, lengths, directions = b
    F_ = MMatrix{3, 4, Float64}(undef)
    torque = MVector{3,Float64}(0,0,0)
    F1 = MVector{3,Float64}(0,0,0)
    for (i, (k, l0, dir, bx)) in enumerate(zip(consts, lengths, directions, bonds))
        rhat, d = norm_and_mag(bx-x)
        # transform bond direction according to bead orientation
        v = orientate_vector(dir, sign(q))
        # torque from diff between rest direction v and actual r
        τ, F = bending_torque_and_force(v, rhat, d, K)
        @. F_[:,i] = -F
        torque += τ
        F1 += F
        F1 += spring_force(rhat, d, l0, k)
    end
    return torque, F1, F_
end


@inline function internal_forces_and_torques!(F, F_, torque, lattice, bead_info, consts)
    K = consts.K
    N = length(lattice)

    @inbounds @fastmath @threads for i in 1:N
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        torque[:,i], F[:,i], F_[:,:,i] = bending_and_spring_forces(lattice.x[i], lattice.q[i], bonds, b, K)
    end

    @threads for i in 1:N
        b = bead_info[i]
        for j in 1:lastindex(b.bonds)
            @. F[:, b.bonds[j]] += F_[:, j, i]
        end
    end
end
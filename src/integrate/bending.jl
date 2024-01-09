"""
    torque_and_force(v, r, K)

Calculate the torque and force due to the bond angle bending between the natural direction v and the actual direction r.
"""
@inline @fastmath function torque_and_force(v, r, K)
    d = norm(r)
    r ./= d
    if all(r .≈ v)
        return zeros(3), zeros(3)
    end

    rdotv = min(1, abs(dot(r,v)))
    # τ = K*sin(θ)*n̂
    torque = K*sqrt(1.0-rdotv^2)*normalize(cross(r, v))
    # F = -K*sin(θ)θ̂/|r|
    force = -K/d*(v - rdotv*r)
    return torque, force
end

"""
    angular_forces!(F, i, b1::Bead, beads, bead_info, dirs, K)::MVector{3,Float64}

Calculate 3D torque and force acting on bead `b1` and its neighbours due to the bond angle bending at `b1`. Updates overall force vectors and returns torque on b1.
"""
function angular_forces!(F, i, b1, beads, bead_info, dirs, K)
    @unpack north, east, south, west = bead_info

    torque = MVector{3,Float64}(0,0,0)
    bonds = [north, east, south, west]
    for (bond,dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            # transform bond direction according to bead orientation
            v = orientate_vector(dir, sign(b1.q))
            @fastmath r = beads[bond].x - b1.x
            # torque from diff between rest direction v and actual r
            τ, F_ = torque_and_force(v, r, K)
            @. F[:,bond] -= F_ 
            @. F[:,i] += F_
            torque += τ
        end
    end
    return torque
end

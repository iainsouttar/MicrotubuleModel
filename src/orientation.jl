

@inline function bond_angle(v, r)
    rhat = normalize(r)
    rdotv = min(1.0,dot(rhat,v))
    thetahat = normalize(v - rdotv*rhat)
    theta = acos(rdotv)
    nhat = cross(thetahat, v)

    return theta, thetahat, nhat
end


@inline function torque_and_force(v, r, K)
    θ, θ_hat, n_hat = bond_angle(v, r)
    if iszero(θ)
        return zeros(3), zeros(3)
    end
    @fastmath torque = K*θ*n_hat
    @fastmath force = -K*θ/norm(r)*θ_hat
    return torque, force
end


function angular_forces!(F, b1, lattice, dirs, K)
    lat, south = b1.lat_nn, b1.long_nn
    north = b1.intra_nn

    torque = MVector{3,Float64}(0,0,0)
    bonds = [north, lat[2], south, lat[1]]
    for (bond,dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            # transform bond direction according to bead orientation
            v = orientate_vector(dir, b1.q)
            @fastmath r = lattice[bond].x - b1.x
            # torque from diff between rest direction v and actual r
            τ, F_ = torque_and_force(v, r, K)
            F[:,bond] += F_
            torque += τ
        end
    end
    return torque
end



function bond_angle(v, r)
    rhat = normalize(r)
    rdotv = LinearAlgebra.dot(rhat,v)
    if isapprox(abs(rdotv), 1.0, atol=1e-4)
        return 0.0, [NaN,NaN,NaN], [NaN,NaN,NaN]
    end
    thetahat = normalize(v - rdotv*rhat)
    theta = acos(rdotv)
    nhat = normalize(cross(rhat, v))

    return theta, thetahat, nhat
end


function torque_and_force(q, v, r, K)
    v_ = orientate_vector(v, q)
    θ, θ_hat, n_hat = bond_angle(v_, r)
    if all(isnan, θ_hat)
        return zeros(3), zeros(3)
    end
    torque = -K*θ*n_hat
    force = -K*θ/norm(r)*θ_hat
    return torque, force
end


function angular_forces!(F, b1, lattice, dirs, K)
    lat, long = b1.lat_nn, b1.long_nn
    intra = b1.intra_nn

    torque = MVector{3,Float64}(0,0,0)
    bonds = [intra,lat[2],long, lat[1]]
    for (bond,dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            τ, F_ = torque_and_force(b1.q, dir, lattice[bond].x - b1.x, K)
            F[:,bond] += F_
            torque += τ
        end
    end
    return torque
end

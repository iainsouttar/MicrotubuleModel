struct BondAngle
    θ::Float64
    ϕ::Float64
end

BondAngle(v::SVector{2,Float64}) = BondAngle(v[1],v[2])

const BondDirec = SVector{3, Float64}

function bond_directions(thetas::Union{AlphaConfirm,BetaConfirm})
    vs = [direc_from_angles(BondAngle(t)) for t in [thetas.north,thetas.east,thetas.south,thetas.west]]
    return SMatrix{3,4, Float64}(hcat(vs...))
end


function quat_from_axisangle(axis::AbstractVector, theta::Real)
    s, c = sincos(theta / 2)
    axis = normalize(axis)
    return Quaternions.Quaternion(c, s*axis[1], s*axis[2], s*axis[3])
end

function to_axis_angle(q::Quaternions.Quaternion)
    r = sqrt(sum(x*x for x in imag_part(q)))
    return BondDirec(imag_part(q)...) ./ r , 2*atan(r,real(q))
end

function direc_from_angles(t::BondAngle)
    x = cos(t.θ)*sin(t.ϕ)
    y = sin(t.θ)*sin(t.ϕ)
    z = cos(t.ϕ)
    return BondDirec(x, y, z)
end

function transform_orientation(v::BondDirec, q::Quaternions.Quaternion)
    v_ = quat(0,v...)
    return conj(q) * v_ * q
end

function orientate_vector(v::BondDirec, q::Quaternions.Quaternion)
    q_prime = transform_orientation(v, q)
    return LinearAlgebra.normalize(BondDirec(imag_part(q_prime)))
end

function bond_angle(v, r)
    rhat = normalize(r)
    #@info norm(rhat), norm(v)
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
    torque = K*θ*n_hat
    force = -K*θ/norm(r)*θ_hat
    return torque, force
end


function angular_forces!(F, b1, lattice, dirs, K)
    lat, long = b1.lat_nn, b1.long_nn
    intra = b1.intra_nn

    torque = MVector{3,Float64}(0,0,0)
    bonds = [intra,lat[1],long, lat[2]]
    for (bond,dir) in zip(bonds,eachcol(dirs))
        if bond != 0
            τ, F_ = torque_and_force(b1.q, dir, lattice[bond].x - b1.x, K)
            F[:,bond] += F_
            torque += τ
        end
    end
    return torque
end

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
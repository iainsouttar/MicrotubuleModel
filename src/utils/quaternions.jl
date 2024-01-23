struct BondAngle
    θ::Float64
    ϕ::Float64
end

BondAngle(v::SVector{2,Float64}) = BondAngle(v[1],v[2])

const BondDirec = SVector{3, Float64}

function direc_from_angles(t::BondAngle)
    s_θ, c_θ = sincos(t.θ)
    s_ϕ, c_ϕ = sincos(t.ϕ)
    return BondDirec(c_θ*s_ϕ, s_θ*s_ϕ, c_ϕ)
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

function transform_orientation(v, q::Quaternions.Quaternion)
    v_ = quat(0, v...)
    return @fastmath conj(q) * v_ * q
end

function orientate_vector(v, q::Quaternions.Quaternion)
    q_prime = transform_orientation(v, q)
    return @fastmath normalize(BondDirec(imag_part(q_prime)))
end
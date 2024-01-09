"""
    iterate!(
        beads::Vector{Bead},
        bead_info::Vector{BeadPars},
        dirs::Dict{Bool,SMatrix{3,4,Float64}},
        conf::Union{PatchConfig, RotationConfig},
        iter::StochEulerPars
    )

Iterate the bead positions and orientations.

Arguments:
- `beads::Vector{Bead}`: mutable position and orientation info
- `bead_info::Vector{BeadPars}`: lattice connections and tubulin type (alpha/beta)
- `dirs::Dict{Bool,SMatrix{3,4,Float64}}`: natural bond angles 
- `conf::Union{PatchConfig, RotationConfig}`: simulation parameters
- `iter`: determines the integration scheme used.

"""
function iterate!(
    beads::Vector{Bead},
    bead_info::Vector{BeadPars},
    #dirs::Dict{Bool,SMatrix{3,4,Float64}},
    conf::Union{PatchConfig,RotationConfig},
    iter::StochEulerPars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    N_tot = lastindex(beads)
    @unpack dt, damp_x, Tk_B = iter
    σ_x = sqrt(dt*2*Tk_B*damp_x)
    F = zeros(Float64, (3, N_tot))
    torque = similar(F)
    
    Δx, Δq = forward(beads, bead_info, F, torque, conf)
    ΔW = σ_x*randn((3,N_tot))
    @inbounds @fastmath @threads for i in 1:N_tot
        beads[i].x += (i>N)*(Δx[i]*dt .+ ΔW[:,i])
        beads[i].q += (i>N)*Δq[i]*dt
    end
end

function iterate!(
    beads::Vector{Bead},
    bead_info::Vector{BeadPars},
    #dirs::Dict{Bool,SMatrix{3,4,Float64}},
    conf::Union{PatchConfig,RotationConfig},
    iter::EulerPars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    F = zeros(Float64, (3, lastindex(beads)))
    F_ = zeros(Float64, (3, 4, lastindex(beads)))
    torque = similar(F)
    
    Δx, Δq = forward(beads, bead_info, F, F_, torque, conf)

    @inbounds @fastmath @threads for i in 1:lastindex(beads)
        beads[i].x += (i>N)*Δx[i]*dt
        beads[i].q += (i>N)*Δq[i]*dt
    end
end

@fastmath @inbounds function iterate!(
    beads::Vector{Bead},
    bead_info::Vector{BeadPars},
    #dirs::Dict{Bool,SMatrix{3,4,Float64}},
    conf::Union{PatchConfig,RotationConfig},
    iter::RK4Pars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = lastindex(beads)
    F = zeros(Float64, (3, lastindex(beads)))
    F_ = zeros(Float64, (3, 4, lastindex(beads)))
    torque = similar(F)
    b = deepcopy(beads)

    N *= 1
    
    k1x, k1q = forward(b, bead_info, F, F_, torque, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k1x[i]*dt/2
        b[i].q = beads[i].q + (i>N)*k1q[i]*dt/2
    end
    k2x, k2q = forward(b, bead_info, F, F_, torque, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k2x[i]*dt/2
        b[i].q = beads[i].q + (i>N)*k2q[i]*dt/2
    end
    k3x, k3q = forward(b, bead_info, F, F_, torque, conf)
    @threads for i in 1:Ntot
        b[i].x = beads[i].x + (i>N)*k3x[i]*dt
        b[i].q = beads[i].q + (i>N)*k3q[i]*dt
    end
    k4x, k4q = forward(b, bead_info, F, F_, torque, conf)

    k_x = @. (k1x + 2*k2x + 2*k3x + k4x)*dt/6
    k_q = @. (k1q + 2*k2q + 2*k3q + k4q)*dt/6
    @threads for i in 1:Ntot
        beads[i].x = beads[i].x + (i>N)*k_x[i]
        beads[i].q = beads[i].q + (i>N)*k_q[i]
    end
end


"""
    forward(
        beads::Vector{Bead},
        bead_info::Vector{BeadPars},
        F::Matrix{Float64},
        torque::Matrix{Float64},
        dirs::Dict{Bool,Union{AlphaConfirm, BetaConfirm}}, 
        conf::Union{PatchConfig, RotationConfig},
    )

Iterate the bead positions and orientations.

Arguments:
- `beads::Vector{Bead}`: mutable position and orientation info
- `bead_info::Vector{BeadPars}`: lattice connections and tubulin type (alpha/beta)
- `F::Matrix{Float64}`: 3D force for each bead to be calculated
- `torque::Matrix{Float64}`: 3D torques for each bead to be calculated
- `dirs::Dict{Bool,SMatrix{3,4,Float64}}`: natural bond angles 
- `conf::Union{PatchConfig, RotationConfig}`: simulation parameters

"""
function forward(
    beads::Vector{Bead},
    bead_info::Vector{BeadPars},
    F::Matrix{Float64},
    F_,
    torque::Matrix{Float64},
    #dirs::Dict{Bool,SMatrix{3,4,Float64}}, 
    conf::Union{PatchConfig, RotationConfig},
)
    Ntot = lastindex(beads)
    x = Vector{BeadPos}(undef, Ntot)
    q = Vector{Quaternions.Quaternion}(undef, Ntot)
    @unpack damp_x, damp_theta, dt = conf.iter_pars

    internal_forces_and_torques!(F, F_, torque, beads, bead_info, conf.spring_consts)

    external_forces!(F, beads, bead_info, conf.external_force)

    @inbounds @fastmath @threads for i in 1:Ntot
        x[i], q[i] = calculcate_deltas(beads[i], F[:,i], torque[:,i], damp_x, damp_theta)
    end
    return x, q
end


@inline function internal_forces_and_torques!(F, F_, torque, beads, bead_info, consts)
    K = consts.K
    N = lastindex(beads)

    @inbounds @fastmath @threads for i in 1:N
        b = bead_info[i]
        bonds = beads[b.bonds]
        torque[:,i], F[:,i], F_[:,:,i] = angular_forces(beads[i], bonds, b.directions, K)
        F[:,i] += linear_forces(beads[i], bonds, bead_info[i])
    end

    @threads for i in 1:N
        b = bead_info[i]
        for j in 1:lastindex(b.bonds)
            @. F[:, b.bonds[j]] += F_[:, j, i]
        end
    end
end

# @inline function internal_forces_and_torques!(F, torque, beads, bead_info, dirs, consts)
#     K = consts.K
#     @inbounds @fastmath @threads for i in 1:lastindex(beads)
#         F[:,i] = linear_forces(beads[i], beads, bead_info[i], consts)
#     end
#     for i in 1:lastindex(beads)
#         torque[:,i] = angular_forces!(F, i, beads[i], beads, bead_info[i], dirs[bead_info[i].α], K)
#     end
# end

"""
    calculcate_deltas(b, F, τ, γ_x, γ_θ)

Calculate the change in the position and orientation for the next timestep.
Returns Δx, Δq.
Ensures the orientation is normalized first

"""
@inline function calculcate_deltas(b, F, τ, γ_x, γ_θ)
    q = sign(b.q)
    q_τ = quat(0, τ...)
    spin = 0.5 * q * q_τ 
    return F * γ_x, spin * γ_θ
end

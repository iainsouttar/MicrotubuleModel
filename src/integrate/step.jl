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
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,RotationConfig},
    iter::StochEulerPars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    Ntot = length(lattice)
    @unpack dt, damp_x, Tk_B = iter
    σ_x = sqrt(dt*2*Tk_B*damp_x)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)

    updateable = (collect(1:Ntot) .> N)
    
    Δx, Δq = forward(lattice, bead_info, F, F_, torque, conf)
    ΔW = σ_x*randn((3,Ntot))
    @inbounds @fastmath @threads for i in 1:Ntot
        lattice.x[i] .+= updateable[i]*(Δx[:,i]*dt .+ ΔW[:,i])
        lattice.q[i] += updateable[i]*Δq[i]*dt
    end
end

function iterate!(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,RotationConfig},
    iter::EulerPars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = length(lattice)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)
    
    Δx, Δq = forward(lattice, bead_info, F, F_, torque, conf)
    updateable = (collect(1:Ntot) .> N)
    #updateable = ones(Bool, Ntot)
    #updateable[1] = false
    update!(lattice, updateable, Δx, Δq, dt)
end

function iterate!(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,RotationConfig},
    iter::RK4Pars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = length(lattice)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)
    lattice_cpy = deepcopy(lattice)

    N *= 1

    updateable = ((1:Ntot) .> N)
    
    k1x, k1q = forward(lattice_cpy, bead_info, F, F_, torque, conf)
    update!(lattice_cpy, lattice, updateable, k1x, k1q, dt/2)

    k2x, k2q = forward(lattice_cpy, bead_info, F, F_, torque, conf)
    update!(lattice_cpy, lattice, updateable, k2x, k2q, dt/2)

    k3x, k3q = forward(lattice_cpy, bead_info, F, F_, torque, conf)
    update!(lattice_cpy, lattice, updateable, k3x, k3q, dt)

    k4x, k4q = forward(lattice_cpy, bead_info, F, F_, torque, conf)
    k_x = @. (k1x + 2*k2x + 2*k3x + k4x)/6
    k_q = @. (k1q + 2*k2q + 2*k3q + k4q)/6
    update!(lattice, updateable, k_x, k_q, dt)
end

@inline function update!(lattice, updateable, Δx, Δq, dt)
    @fastmath @inbounds @threads for i in 1:length(lattice)
        lattice.x[i] .+= updateable[i]*Δx[:,i]*dt
        lattice.q[i] += updateable[i]*Δq[i]*dt
    end
end

@inline function update!(lattice_cpy, lattice, updateable, Δx, Δq, dt)
    @fastmath @inbounds @threads for i in 1:length(lattice)
        lattice_cpy.x[i] .= lattice.x[i] .+ updateable[i]*Δx[:,i]*dt
        lattice_cpy.q[i] = lattice.q[i] + updateable[i]*Δq[i]*dt
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
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    F::Matrix{Float64},
    F_,
    torque::Matrix{Float64},
    conf::Union{PatchConfig, RotationConfig},
)
    Ntot = length(lattice)
    q = Vector{Quaternions.Quaternion}(undef, Ntot)
    @unpack damp_x, damp_theta, dt = conf.iter_pars

    internal_forces_and_torques!(F, F_, torque, lattice, bead_info, conf.spring_consts)

    external_forces!(F, lattice, bead_info, conf.external_force)

    @inbounds @fastmath @threads for i in 1:Ntot
        q[i] = calculcate_deltas(lattice.q[i], torque[:,i], damp_theta)
    end
    return F.*damp_x, q
end


@inline function internal_forces_and_torques!(F, F_, torque, lattice, bead_info, consts)
    K = consts.K
    N = length(lattice)

    # this loop could possibly be much faster if distributed over many threads
    @inbounds @fastmath @threads for i in 1:N
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        #F[:, i] .= 0.0
        torque[:,i], F[:,i], F_[:,:,i] = angular_forces(lattice.x[i], lattice.q[i], bonds, b.directions, K)
        F[:,i] += linear_forces(lattice.x[i], bonds, bead_info[i])
    end

    @threads for i in 1:N
        b = bead_info[i]
        for j in 1:lastindex(b.bonds)
            @. F[:, b.bonds[j]] += F_[:, j, i]
        end
    end
end

"""
    calculcate_deltas(b, F, τ, γ_x, γ_θ)

Calculate the change in the position and orientation for the next timestep.
Returns Δx, Δq.
Ensures the orientation is normalized first

"""
@inline function calculcate_deltas(q, τ, γ_θ)
    q_τ = quat(0, τ...)
    spin = 0.5 * sign(q) * q_τ 
    return spin * γ_θ
end

"""
    iterate!(
        beads::Lattice,
        bead_info::Vector{BeadPars},
        conf::Union{PatchConfig, MicrotubuleConfig},
        iter::StochEulerPars
    )

Iterate the bead positions and orientations.

Arguments:
- `beads::Lattice`: mutable position and orientation of the beads
- `bead_info::Vector{BeadPars}`: lattice connections and tubulin type (alpha/beta)
- `dirs::Dict{Bool,SMatrix{3,4,Float64}}`: natural bond angles 
- `conf::Union{PatchConfig, MicrotubuleConfig}`: simulation parameters
- `iter`: determines the integration scheme used.

"""
function iterate!(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,MicrotubuleConfig},
    iter::StochEulerPars,
    flag::Int64
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    Ntot = length(lattice.x)
    @unpack dt, damp_x, Tk_B = iter
    σ_x = sqrt(dt*2*Tk_B*damp_x)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)

    updateable = (collect(1:Ntot) .> N)
    #updateable = ones(Bool, Ntot)
    #updateable[1] = false
    
    Δx, Δq = forward(lattice, bead_info, F, F_, torque, conf,flag)
    ΔW = σ_x*randn((3,Ntot))
    @inbounds @fastmath @threads for i in 1:Ntot
        lattice.x[i] .+= updateable[i]*(Δx[:,i]*dt .+ ΔW[:,i])
        lattice.q[i] += updateable[i]*Δq[i]*dt
        if lattice.kinesin[i,1] == false
            random_num = rand(Exponential(1/lattice.rates[i]))
            if (random_num <= dt) & (i>2*N)
                make_kinesin_like!(lattice, bead_info, conf.lattice.N, conf.lattice.S, i)
                #print("changed ", i)
            end
        end

    end
end

function iterate!(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,MicrotubuleConfig},
    iter::EulerPars,
    flag::Int64
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = length(lattice.x)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)
    
    Δx, Δq = forward(lattice, bead_info, F, F_, torque, conf,flag)


    #this choiceseems to actually make a difference
    #print(F[:,1], torque[:,1], Δx[:,1], Δq[:,1], N, "\n")
    updateable = (collect(1:Ntot) .> N)
    #updateable = ones(Bool, Ntot)
    #updateable[1] = false
    update!(lattice, updateable, conf, bead_info, Δx, Δq, dt)
end

function iterate!(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    conf::Union{PatchConfig,MicrotubuleConfig},
    iter::RK4Pars
)
    N = isa(conf.lattice,LatticePatchPars) ? conf.lattice.N_lat : conf.lattice.N
    dt = iter.dt
    Ntot = length(lattice.x)
    F = zeros(Float64, (3, Ntot))
    F_ = zeros(Float64, (3, 4, Ntot))
    torque = similar(F)
    lattice_cpy = deepcopy(lattice)

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


@inline function update!(lattice, updateable, conf, bead_info, Δx, Δq, dt)
    @fastmath @inbounds @threads for i in 1:length(lattice.x)
        #print(updateable[i])
        lattice.x[i] .+= updateable[i]*Δx[:,i]*dt
        lattice.q[i] += updateable[i]*Δq[i]*dt
        if lattice.kinesin[i,1] == false
            random_num = rand(Exponential(1/lattice.rates[i]))
            if random_num <= dt
                make_kinesin_like!(lattice, bead_info, conf.lattice.N, conf.lattice.S, i)
                #print("changed ", i)
            end
        end
        
    end
end

function simulate_event(lambda, interval_length)
    # Generate a random number from exponential distribution
    random_number = rand(Exponential(1/lambda))
    
    # Check if the random number is within the interval
    return random_number <= interval_length
end

@inline function update!(lattice_cpy, lattice, updateable, Δx, Δq, dt)
    @fastmath @inbounds @threads for i in 1:length(lattice.x)
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
        conf::Union{PatchConfig, MicrotubuleConfig},
    )

Iterate the bead positions and orientations.

Arguments:
- `beads::Vector{Bead}`: mutable position and orientation info
- `bead_info::Vector{BeadPars}`: lattice connections and tubulin type (alpha/beta)
- `F::Matrix{Float64}`: 3D force for each bead to be calculated
- `torque::Matrix{Float64}`: 3D torques for each bead to be calculated
- `dirs::Dict{Bool,SMatrix{3,4,Float64}}`: natural bond angles 
- `conf::Union{PatchConfig, MicrotubuleConfig}`: simulation parameters

"""
function forward(
    lattice::Lattice,
    bead_info::Vector{BeadPars},
    F::Matrix{Float64},
    F_,
    torque::Matrix{Float64},
    conf::Union{PatchConfig, MicrotubuleConfig},
    flag::Int64
)
    Ntot = length(lattice.x)
    q = Vector{Quaternions.Quaternion}(undef, Ntot)
    @unpack damp_x, damp_theta, dt = conf.iter_pars

    internal_forces_and_torques!(F, F_, torque, lattice, bead_info,flag, conf.lattice.N, conf.lattice.S)
    external_forces!(F, lattice, bead_info, conf.external_force)
    
    @inbounds @fastmath @threads for i in 1:Ntot
        q[i] = quat_delta(lattice.q[i], torque[:,i], damp_theta)
    end
    #q[1] = quat_from_axisangle([0,0,1],-π/2)
    #F[:,1] = MVector{3,Float64}(0,0,0)



    return F.*damp_x, q
end

"""
    quat_delta(q, τ, γ_θ)

Calculate the change in the position and orientation for the next timestep.
"""
quat_delta(q, τ, γ_θ) = γ_θ * 0.5 * sign(q) * quat(0, τ...)

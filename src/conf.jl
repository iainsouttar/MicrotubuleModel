@option "iterate" struct IteratePars
    dt::Float64
    mass::Float64
    moment_inertia::Float64
    damp_x::Float64
    damp_theta::Float64
end

@option "lattice" struct LatticePars
    S::Int
    N::Int
    num_rings::Int
    a::Float64
    dx::Float64
end

@option "patch" struct LatticePatchPars
    S::Int
    N::Int
    N_lat::Int
    N_long::Int
    a::Float64
    dx::Float64
end

@option "spring consts" struct SpringConst
    K::Float64
    k_long::Float64
    k_lat::Float64
    k_in::Float64
    k_in_kin::Float64
    l0_long::Float64
    l0_lat::Float64
    l0_in::Float64
    l0_in_kin::Float64
end

@option "alpha struct" struct AlphaConfirm
    north::SVector{2,Float64} = [π/2, -0.2]
    east::SVector{2,Float64} = [π/13, π/2-0.1819]
    south::SVector{2,Float64} = [π/2, π]
    west::SVector{2,Float64} = [π-π/13,π/2+0.1819]
end

@option "beta struct" struct BetaConfirm
    north::SVector{2,Float64} = [π/2, 0]
    east::SVector{2,Float64} = [π/13, π/2-0.1819]
    south::SVector{2,Float64} = [π/2, π+0.2]
    west::SVector{2,Float64} = [π-π/13,π/2+0.1819]
end

"""
Full lattice and simulation parameters
"""
@option "rotation" struct RotationConfig
    lattice::LatticePars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::IteratePars
    spring_consts::SpringConst
end

"""
Parameters for a patch of a lattice e.g. 5x5 grid 
"""
@option "patch conf" struct PatchConfig
    lattice::LatticePatchPars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::IteratePars
    spring_consts::SpringConst
end
@option "euler" struct EulerPars
    dt::Float64
    damp_x::Float64
    damp_theta::Float64
end

@option "stoch_euler" struct StochEulerPars
    dt::Float64
    Tk_B::Float64
    damp_x::Float64
    damp_theta::Float64
end

@option "rk4" struct RK4Pars
    dt::Float64
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

@option "istropic_force" struct IsotropicForce
    F::Float64 = 1.0
end

@option "youngs_modulus" struct YoungsModulusTest
    F::Float64 = 1.0
    N::Int = 13
end

@option "none" struct NoExternalForce end

"""
Full lattice and simulation parameters
"""
@option "rotation" struct RotationConfig
    lattice::LatticePars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::Union{RK4Pars, EulerPars, StochEulerPars}
    spring_consts::SpringConst
    external_force::Union{NoExternalForce, IsotropicForce, YoungsModulusTest}
end

"""
Parameters for a patch of a lattice e.g. 5x5 grid 
"""
@option "patch conf" struct PatchConfig
    lattice::LatticePatchPars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::Union{RK4Pars, EulerPars, StochEulerPars}
    spring_consts::SpringConst
    external_force::Union{NoExternalForce, IsotropicForce, YoungsModulusTest}
end

"""
Full lattice and simulation parameters
"""
@option "YM conf" struct YoungsModulus
    lattice::LatticePars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::Union{RK4Pars, EulerPars, StochEulerPars}
    spring_consts::SpringConst
    external_force::Union{NoExternalForce, YoungsModulusTest}
end
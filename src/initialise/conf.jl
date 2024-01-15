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

####################################################

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

###########################################################

# External forces for various mechanical tests

@option "istropic_force" struct IsotropicForce
    F::Float64 = 1.0
end

@option "youngs_modulus" struct YoungsModulusTest
    F::Float64 = 1.0
    N::Int = 13
end

@option "youngs_modulus_pf" struct YoungsModulusPF
    F::Float64 = 1.0
    N::Int = 13
end

@option "bending_stiffness" struct BendingStiffnessTest
    F::Float64 = 1.0
    N::Int = 13
end

@option "none" struct NoExternalForce end


########################################################

function set_bond_angles(conf)
    @unpack S, N, dx, a = conf.lattice
    α, β = calc_natural_angles(S, N, dx, a)
    conf = @set conf.alpha = α
    conf = @set conf.beta = β
    return conf
end

"""
Full lattice and simulation parameters
"""
@option "rotation" struct MicrotubuleConfig
    lattice::LatticePars
    alpha::AlphaConfirm = AlphaConfirm()
    beta::BetaConfirm = BetaConfirm()
    iter_pars::Union{RK4Pars, EulerPars, StochEulerPars}
    spring_consts::SpringConst
    external_force::Union{
        NoExternalForce, 
        IsotropicForce, 
        YoungsModulusTest,
        YoungsModulusPF,
        BendingStiffnessTest}
end

"""
    initialise(conf)::Tuple{Lattice, Vector{BeadPars}}

Initialise the lattice and construct the bond directions
"""
function initialise(conf::MicrotubuleConfig)::Tuple{Lattice, Vector{BeadPars}}
    @unpack num_rings, S, N, a, dx = conf.lattice

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    beads, structure = create_lattice(dirs, conf.spring_consts, num_rings, a, dx; S=S, N=N)

    return beads, structure

end

"""
    initialise_dimer(conf)::Tuple{Lattice, Vector{BeadPars}}

Initialise the dimer as a lattice and construct the bond directions
"""
function initialise_dimer(conf::MicrotubuleConfig)
    @unpack S, N, a, dx = conf.lattice

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    beads, structure = create_dimer(dirs, conf.spring_consts, a, dx; S=S, N=N)

    return beads, structure

end


"""
    initialise_PF(conf)::Tuple{Lattice, Vector{BeadPars}}

Initialise the protofilament as a lattice and construct the bond directions
"""
function initialise_PF(conf::MicrotubuleConfig)
    @unpack S, N, a, dx, num_rings = conf.lattice

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    beads, structure = create_PF(dirs, conf.spring_consts, num_rings, a, dx; S=S, N=N)

    return beads, structure

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
    external_force::Union{
        NoExternalForce, 
        IsotropicForce, 
        YoungsModulusTest, 
        BendingStiffnessTest}
end

function initialise(conf::PatchConfig)
    @unpack N_lat, N_long, S, N, a, dx = conf.lattice

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    beads, structure = create_patch(dirs, conf.spring_consts, N_lat, N_long, a, dx; S=S, N=N)

    return beads, structure

end

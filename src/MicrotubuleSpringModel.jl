module MicrotubuleSpringModel

using StaticArrays
using LinearAlgebra
using LinearAlgebra: norm, normalize
using Parameters
using Base.Threads
using Distributions
using Random
using Quaternions
using Configurations: @option
using Setfield: @set
using LoopVectorization: @tturbo

using GLMakie
using ColorSchemes

export 
    BeadPos,
    BeadAngle,
    Bead,
    create_lattice,
    bond_directions,
    neighbours,
    set_bond_angles,

    total_energy,
    LinSpringConst,

    quat_from_axisangle,
    orientate_vector,

    BondAngle,
    BondDirec,
    torque_and_force,

    iterate!,
    iterateSDE!,

    RotationConfig,
    PatchConfig,

    surface_area,
    youngs_modulus,
    microtubule_length,

    plot,
    colorschemes


BeadPos = MVector{3, Float64}
BeadAngle = MVector{3, Float64}

# three vector for position and for orientation angles
# alpha true for alpha monomer, false for beta monomer 
mutable struct Bead
    x::BeadPos
    q::Quaternions.Quaternion
    kinesin::Bool
end

# three vector for position and for orientation angles
# alpha true for alpha monomer, false for beta monomer 
struct BeadPars
    α::Bool
    north::Int
    east::Int
    south::Int
    west::Int
end

include("conf.jl")
include("transformations.jl")
include("lattice.jl")
include("integrate/forces.jl")
include("integrate/energy.jl")
include("integrate/step.jl")
include("integrate/orientation.jl")
include("visuals.jl")
include("utils.jl")

"""
    initialise(conf)::Tuple{Vector{Bead}, Dict{Bool,BondDirec}}

    Initialise the lattice and construct the bond directions
"""
function initialise(conf::RotationConfig)::Tuple{Vector{Bead}, Vector{BeadPars}, Dict{Bool,SMatrix{3,4, Float64}}} 
    @unpack num_rings, S, N, a, dx = conf.lattice
    beads, structure = create_lattice(num_rings, a, dx; S=S, N=N)

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    return beads, structure, dirs

end

function initialise(conf::PatchConfig)
    @unpack N_lat, N_long, S, N, a, dx = conf.lattice
    beads, structure = create_patch(N_lat, N_long, a, dx; S=S, N=N)

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    return beads, structure, dirs

end


function initialise_dimer(conf::RotationConfig)
    @unpack S, N, a, dx = conf.lattice
    beads, structure = create_dimer(a, dx; S=S, N=N)

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    return beads, structure, dirs

end

end

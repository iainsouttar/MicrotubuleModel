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

    plot,
    colorschemes


BeadPos = MVector{3, Float64}
BeadAngle = MVector{3, Float64}

# three vector for position and for orientation angles
# alpha true for alpha monomer, false for beta monomer 
mutable struct Bead
    x::BeadPos
    q::Quaternions.Quaternion
    α::Bool
    kinesin::Bool
    lat_nn::Tuple{Int,Int}
    long_nn::Int
    intra_nn::Int
end

include("conf.jl")
include("lattice.jl")
include("forces.jl")
include("transformations.jl")
include("orientation.jl")
include("sde.jl")
include("visuals.jl")

"""
    initialise(conf)::Tuple{Vector{Bead}, Dict{Bool,BondDirec}}

    Initialise the lattice and construct the bond directions
"""
function initialise(conf::RotationConfig)::Tuple{Vector{Bead}, Dict{Bool,SMatrix{3,4, Float64}}} 
    @unpack num_rings, S, N, a, dx = conf.lattice
    lattice = create_lattice(num_rings, a, dx; S=S, N=N)

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    return lattice, dirs

end

function initialise(conf::PatchConfig)
    @unpack N_lat, N_long, S, N, a, dx = conf.lattice
    lattice = create_patch(N_lat, N_long, a, dx; S=S, N=N)

    dirs = Dict(
        true => bond_directions(conf.alpha),
        false => bond_directions(conf.beta)
    )

    return lattice, dirs

end

end

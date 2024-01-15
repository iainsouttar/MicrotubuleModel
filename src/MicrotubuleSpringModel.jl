module MicrotubuleSpringModel

using StaticArrays
using LinearAlgebra
using LinearAlgebra: norm, normalize
using Parameters
using Base.Threads
using Distributions
using Random
using Quaternions
using Configurations: @option, to_toml
using Setfield: @set
using DelimitedFiles
using CSV
using DataFrames
using ProgressMeter: @showprogress

using GLMakie
using ColorSchemes

export 
    BeadPos,

    Lattice,
    BeadPars,
    SpringConst,

    create_lattice,
    initialise,
    initialise_dimer,
    initialise_PF,
    set_bond_angles,

    total_energy,

    iterate!,

    MicrotubuleConfig,
    PatchConfig,

    surface_area,
    youngs_modulus,
    microtubule_length,
    deflection_end,
    stiffness,

    save_to_csv,
    save_params,
    load_from_csv,

    plot,
    plot_flat!,
    plot_individual!,
    colorschemes,

    burnin


BeadPos = MVector{3, Float64}

mutable struct Lattice{N}
    x::MVector{N, BeadPos}
    q::Vector{Quaternions.Quaternion}
    kinesin::MVector{N, Bool}
end

function Lattice(x, q, kinesin)
    N = length(kinesin)
    return Lattice{N}(x, q, kinesin)
end

Base.size(l::Lattice) = size(l.kinesin)

Base.length(l::Lattice) = length(l.kinesin)


struct BeadPars
    α::Bool
    bonds::Vector{Int}
    directions::Vector{SVector{3,Float64}}
    consts::Vector{Float64}
    lengths::Vector{Float64}
end

include("utils/quaternions.jl")

include("initialise/lattice.jl")
include("initialise/neighbours.jl")
include("initialise/angles.jl")
include("initialise/conf.jl")

include("integrate/springs.jl")
include("integrate/bending.jl")
include("integrate/external_forces.jl")
include("integrate/energy.jl")
include("integrate/step.jl")

include("simulate.jl")

include("utils/visuals.jl")
include("utils/measurements.jl")
include("utils/io.jl")

end

module MicrotubuleSpringModel

using StaticArrays
using LinearAlgebra
using Parameters
using Base.Threads
using Distributions
using Random

using GLMakie
using ColorSchemes

export 
    BeadPos,
    BeadAngle,
    Bead,
    create_lattice,
    neighbours,

    total_energy,
    LinSpringConst,
    iterate!,
    iterateSDE!,

    plot


BeadPos = MVector{3, Float64}
BeadAngle = MVector{3, Float64}

# three vector for position and for orientation angles
# alpha true for alpha monomer, false for beta monomer 
mutable struct Bead
    x::BeadPos
    θ::BeadAngle
    α::Bool
    # lat_nn::Tuple{Ref,Ref}
    # long_nn::Ref
    # intra_nn::Ref
    lat_nn::Tuple{Int,Int}
    long_nn::Int
    intra_nn::Int
end

include("lattice.jl")
include("forces.jl")
include("sde.jl")
include("visuals.jl")

end

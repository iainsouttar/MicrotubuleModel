#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using StaticArrays
using Logging
using LinearAlgebra
using LinearAlgebra: normalize, norm
using Distributions
using Parameters
using Colors
using GLMakie

using Quaternions
using Quaternions: Quaternion

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/bending_stiffness.toml")

conf = set_bond_angles(conf)

beads, bead_info = MicrotubuleSpringModel.initialise(conf)

@benchmark iterate!($beads, $bead_info, $conf, $conf.iter_pars) samples=500

#######################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stochastic.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

@benchmark iterate!($beads, $bead_info, $dirs, $conf, $conf.iter_pars)

###############################################

b1 = beads[15]
@unpack north, α = bead_info[15]

v = MicrotubuleSpringModel.orientate_vector(dirs[α][:,1], b1.q)
@fastmath r = beads[north].x - b1.x

@benchmark torque_and_force($v, $r, 1.0)
#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
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

conf = from_toml(MicrotubuleConfig, "config/bending_stiffness.toml")

conf = set_bond_angles(conf)

beads, bead_info = initialise(conf)

@benchmark iterate!($beads, $bead_info, $conf, $conf.iter_pars)

#######################################################

conf = from_toml(MicrotubuleConfig, "config/equilibrium.toml")

conf = set_bond_angles(conf)

lattice, bead_info = initialise(conf)

@benchmark iterate!($lattice, $bead_info, $conf, $conf.iter_pars)



conf = from_toml(MicrotubuleConfig, "config/equilibrium.toml")
conf = set_bond_angles(conf)

beads, bead_info = MicrotubuleSpringModel.initialise(conf)

Nt = Int(1e6)

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, conf, conf.iter_pars)
end

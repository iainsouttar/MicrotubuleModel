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

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stochastic.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

Nt = 10000
step = 1
time = 0:step:Nt
E = zeros((6,Nt÷step+1))

E[:,1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    if i % step == 0
        E[:,i÷step+1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
    end
end

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)
scene
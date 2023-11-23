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

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")

conf = set_bond_angles(conf)

lattice, dirs = MicrotubuleSpringModel.initialise(conf)

# for (i,b) in enumerate(lattice)
#     if i % 13 ∈ (4,5,6,7,8)
#         b.kinesin = true
#     end
# end

@showprogress for i in 1:10000
    iterate!(lattice, conf, dirs)
end

GLMakie.closeall()
scene = plot(lattice, dirs)
scene



norms = [norm(v) for v in eachcol(dirs[false])]


k_ab = 3.5
k_bb = 1.0/0.4
k_a = 3.6
k_b = 3.8
K_ab = 1/(1/(2*k_a)+1/(2*k_b)+1/k_ab)

K_bb = 1/(1/k_b+1/k_bb)
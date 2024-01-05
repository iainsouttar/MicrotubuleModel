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

Nt = 100000
s = zeros(Float64,(4,Nt))
E = zeros(Nt)

@showprogress for i in 1:Nt
    iterate!(lattice, conf, dirs)
    E[i] = total_energy(lattice, dirs, conf.spring_consts)
end

nhat = zeros(Float64, (3,4,lastindex(lattice)))
nhat[1,:,:] .= 1

for (i,b) in enumerate(lattice)
    lat, long = b.lat_nn, b.long_nn
    intra = b.intra_nn
    bonds = [intra,lat[2],long, lat[1]]
    for (j,(bond,dir)) in enumerate(zip(bonds,eachcol(dirs[b.α])))
        if bond != 0
            # transform bond direction according to bead orientation
            v = MicrotubuleSpringModel.orientate_vector(dir, b.q)
            r = lattice[bond].x - b.x
            # torque from diff between rest direction v and actual r
            #θ, θ_hat, nhat[:,j, i] = MicrotubuleSpringModel.bond_angle(v, r)
            τ, F_ = torque_and_force(v, r, 1.0)
            nhat[:,j, i] .= 2*τ
        end
    end
end

GLMakie.activate!()
GLMakie.closeall()
scene = plot(lattice, dirs)


xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([2*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)

scene

F = zeros(Float64, (3, lastindex(lattice)))
torque = similar(F)

MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(lattice, dirs)

xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([10*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)


scene
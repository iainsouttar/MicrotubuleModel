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

Nt = 20_000
s = zeros(Float64,(4,Nt))

@showprogress for i in 1:Nt
    iterate!(lattice, conf, dirs)
    for b in lattice 
        lat, long = b.lat_nn, b.long_nn
        intra = b.intra_nn
        bonds = [intra,lat[2],long, lat[1]]
        for (j,(bond,dir)) in enumerate(zip(bonds,eachcol(dirs[b.α])))
            if bond != 0
                v_ = MicrotubuleSpringModel.orientate_vector(dir, b.q)
                θ, θ_hat, n_hat = MicrotubuleSpringModel.bond_angle(v_, lattice[bond].x - b.x)
                s[j,i] += abs(θ)
            end

        end
    end

    s[:,i] ./= lastindex(lattice)
end

GLMakie.activate!()
GLMakie.closeall()
scene = plot(lattice, dirs)
idx = 7*13+5
lat = lattice[idx].lat_nn
long = lattice[idx].long_nn
intra = lattice[idx].intra_nn
s = Scene(scene, camera=scene.camera)
lat_plots = [mesh!(s, Sphere(Point3f(lattice[b].x), a/2), color=:black) for b in lat if b!=0]
long_plots = long==0 ? Nothing : mesh!(s, Sphere(Point3f(lattice[long].x), a/2), color=:yellow)
long_plots = intra==0 ? Nothing : mesh!(s, Sphere(Point3f(lattice[intra].x), a/2), color=:purple)
mesh!(s, Sphere(Point3f(lattice[idx].x), a/2), color=:pink)
scene

F = zeros(Float64, (3, lastindex(lattice)))
torque = similar(F)

MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

GLMakie.closeall()
scene = plot(lattice, dirs)

xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([5*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)


scene



norms = [norm(v) for v in eachcol(dirs[false])]


k_ab = 3.5
k_bb = 1.0/0.4
k_a = 3.6
k_b = 3.8
K_ab = 1/(1/(2*k_a)+1/(2*k_b)+1/k_ab)

K_bb = 1/(1/k_b+1/k_bb)


v = BondDirec(1,0,0)
r = BondDirec(0,1,0)

θ, θ_hat, n_hat = MicrotubuleSpringModel.bond_angle(v,r)

v = MicrotubuleSpringModel.direc_from_angles(BondAngle(0,π/2))
r = MicrotubuleSpringModel.direc_from_angles(BondAngle(π/4,π/2))

θ, θ_hat, n_hat = MicrotubuleSpringModel.bond_angle(v,r)

s = [0.0,0.0,0.0,0.0]
for b in lattice 
    lat, long = b.lat_nn, b.long_nn
    intra = b.intra_nn
    bonds = [intra,lat[2],long, lat[1]]
    for (i,(bond,dir)) in enumerate(zip(bonds,eachcol(dirs[b.α])))
        if bond != 0
            v_ = MicrotubuleSpringModel.orientate_vector(dir, b.q)
            θ, θ_hat, n_hat = MicrotubuleSpringModel.bond_angle(v_, lattice[bond].x - b.x)
            s[i] += θ
        end

    end
end

s /= lastindex(lattice)


CairoMakie.activate!()
f = Figure()
ax = Axis(f[1,1])
series!(ax,collect(1:Nt),s)
xlims!(0,Nt)
f
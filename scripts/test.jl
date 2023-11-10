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
using Distributions
using Parameters
using Colors
using GLMakie

struct LatticePars
    S::Int
    N::Int
    num_rings::Int
    a::Float64
    δx::Float64
end

rings = 20
a = 4.05

params = LatticePars(3, 13, rings, a, 5.13)
consts = LinSpringConst(1.0, 1.0, 1.2, 1.5, a, 1.2*a, 0.9*a, 1.3*a)

function main(
    consts::LinSpringConst,
    params::LatticePars,
    steps::Int,
    dt::Real;
    kBT=0.005
)
    @unpack num_rings, S, N, a, δx = params
    lattice = create_lattice(num_rings, a, δx; S=S, N=N)

    for (i,b) in enumerate(lattice)
        if i % 13 ∈ (4,5,6,7,8)
            b.kinesin = true
        end
    end

    energy = zeros(Float64, steps÷50)
    dW = zeros((3,length(lattice)))
    noise = Normal(0.0, sqrt(2*kBT*dt))

    @showprogress for i in 1:steps
        iterateSDE!(lattice, consts, noise, dW, dt)
        if i % 50 == 0
            energy[i÷50] = total_energy(lattice, consts)
        end
    end

    return lattice, energy
end

long_angle_0 = BeadAngle(0.1, 0.0, 0.0)
lat_angle_0 = BeadAngle(0.0, 0.2, 0.4)


lattice, energy = main(consts, params, 10000, 0.01)

GLMakie.activate!(ssao=true)
GLMakie.closeall() # close any open screen

scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
plot!(scene, lattice)
scene

save("kinesin-buckling-2.png", scene, px_per_unit = 2)

@unpack num_rings, S, N, a, δx = params
lattice = create_lattice(num_rings, a, δx; S=S, N=N)

f = Figure()
ax = Axis(f[1,1])
lines!(ax, collect(1:50:steps), energy)
f


using GLMakie

GLMakie.activate!()


lattice = create_lattice(7, a, δx; S=3, N=13)
scene = plot_lattice(lattice)

idx = 3*13+7
lat = lattice[idx].lat_nn
long = lattice[idx].long_nn
intra = lattice[idx].intra_nn

s = Scene(scene, camera=scene.camera)
lat_plots = [mesh!(s, Sphere(Point3f(lattice[b].x), a/2), color=:black) for b in lat if b!=0]
long_plots = long==0 ? Nothing : mesh!(s, Sphere(Point3f(lattice[long].x), a/2), color=:yellow)
long_plots = intra==0 ? Nothing : mesh!(s, Sphere(Point3f(lattice[intra].x), a/2), color=:purple)
mesh!(s, Sphere(Point3f(lattice[idx].x), a/2), color=:pink)

scene

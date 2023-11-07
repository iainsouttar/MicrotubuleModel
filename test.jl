#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using StaticArrays
using Logging
using LinearAlgebra

S = 3
N = 13
a = 4.05
δx = 5.13
num_rings = 10

long_angle_0 = BeadAngle(0.1, 0.0, 0.0)
lat_angle_0 = BeadAngle(0.0, 0.2, 0.4)


consts = LinSpringConst(1.0, 1.0, 1.2, a, 1.2*a, 0.9*a)
lattice = create_lattice(num_rings, a, δx; S=3, N=13)


steps = 10000
energy = zeros(Float64, steps÷50)

@showprogress for i in 1:steps
    iterate!(lattice, consts, 0.05)
    #if i % 50 == 0
        #energy[i÷50] = sum(bead_energy(b, lattice, consts) for (idx,b) in enumerate(lattice))
    #end
end



scene = plot_lattice(lattice)
scene



f = Figure()
ax = Axis(f[1,1])
lines!(ax, collect(1:50:steps), energy)
f


mutable struct Tst3
    x::Int 
    a::Union{Tst3,Nothing}
end

t = Tst3(5,nothing)

a = Tst3(4,nothing)

t.a = a

a.x = 3

t
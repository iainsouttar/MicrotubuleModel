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


energy = zeros(Float64, steps÷50)
dW = zeros((3,length(lattice)))
noise = Normal(0.0, sqrt(2*kBT*dt))

@showprogress for i in 1:steps
    iterateSDE!(lattice, consts, noise, dW, dt)
    if i % 50 == 0
        energy[i÷50] = total_energy(lattice, consts)
    end
end

f = Figure()
ax = Axis(f[1,1])
lines!(ax, collect(1:50:steps), energy)
f

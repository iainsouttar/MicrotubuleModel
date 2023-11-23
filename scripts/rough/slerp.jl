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

using Quaternions
using Quaternions: Quaternion




GLMakie.rotate!(arrs, GLMakie.Quaternion(norm))

scene

slerps = [slerp(dir, norm, t) for t in 0.1:0.2:0.9]

for (i,sl) in enumerate(slerps)
    arrs = arrows!(s, [Point3f(0),Point3f(0),Point3f(0)], Vec3f.(vecs), linewidth=0.1, color=colors[i])
    #arrs = arrows!(s, [Point3f(0)], [Vec3f(1,0,0)], linewidth=0.1, color=colors[i])
end

arrs = arrows!(scene, [Point3f(0)], [Vec3f(-1,0,0)], linewidth=0.1, color=:black)

scene
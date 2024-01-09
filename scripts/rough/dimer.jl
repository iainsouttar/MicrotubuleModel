#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using Quaternions
using Quaternions: Quaternion



conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")
conf = set_bond_angles(conf)
beads, bead_info, dirs = MicrotubuleSpringModel.initialise_dimer(conf)

Nt = 100
@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
end

F = zeros(Float64, (3, lastindex(beads)))
torque = similar(F)
MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, beads, bead_info, dirs, conf.spring_consts)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info, dirs)
xs = [Point3f(b.x) for b in beads]
vs = Vector{Vec{3, Float32}}([2*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)
scene

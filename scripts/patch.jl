if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


conf = from_toml(PatchConfig, "config/patch.toml")
conf = MicrotubuleSpringModel.set_bond_angles(conf)

lattice, dirs = MicrotubuleSpringModel.initialise(conf)

@showprogress for i in 1:10000
    iterate!(lattice, conf, dirs)
end

F = zeros(Float64, (3, lastindex(lattice)))
torque = similar(F)

MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

GLMakie.closeall()
scene = plot(lattice, dirs)

xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([2.5*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)

scene
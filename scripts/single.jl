if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end



q = quat_from_axisangle([0,0,1],0)

bead = Bead(
    BeadPos(0,0,0), 
    q,
    false, 
    false,
    (0,0), 0, 0
)

dt = 0.0004
torque = MVector{3, Float64}(0,1,0)
F = MVector{3, Float64}(0,0,0)
v = BondDirec(1,0,0)

directions = []

for i in 1:5000
    MicrotubuleSpringModel.step!(bead, F, torque, dt, dt)
    if i % 200 == 0
        v_ = orientate_vector(v, bead.q)
        push!(directions, v_)
    end
end

GLMakie.activate!()
GLMakie.closeall()
scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
arrows!(scene, [Point3f(bead.x)], Vector{Vec{3, Float32}}(directions), linewidth=0.05, color=1:length(directions), arrowsize=0.1)
arrows!(scene, [Point3f(bead.x)], Vector{Vec{3, Float32}}([torque]), linewidth=0.2, color=:white, arrowsize=0.3)
mesh!(scene, Sphere(Point3f(bead.x), 0.2), color=MicrotubuleSpringModel.COLORS[(bead.α, bead.kinesin)], shininess=32.0)
GLMakie.scale!(scene, 0.05, 0.05, 0.05)
center!(scene)
return scene





scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)

for b in lattice
    v = MicrotubuleSpringModel.transform_orientation(BondDirec(0,1,0),b.q)
    v_ = 5*normalize([imag_part(v)...])
    arrows!(scene, [Point3f(b.x)], Vector{Vec{3, Float32}}([v_]), linewidth=0.2, color=:white, arrowsize=0.3)
    v = MicrotubuleSpringModel.transform_orientation(BondDirec(1,0,0),b.q)
    v_ = 5*normalize([imag_part(v)...])
    arrows!(scene, [Point3f(b.x)], Vector{Vec{3, Float32}}([v_]), linewidth=0.2, color=MicrotubuleSpringModel.NATURE.colors[4], arrowsize=0.3)
end

for b in lattice
    mesh!(scene, Sphere(Point3f(b.x), 1.25), color=MicrotubuleSpringModel.COLORS[(b.α, b.kinesin)], shininess=32.0)
end

GLMakie.scale!(scene, 0.05, 0.05, 0.05)
center!(scene)
return scene
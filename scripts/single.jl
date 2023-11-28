if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


v = BondDirec(1,0,0)
dt = 0.1

q = MicrotubuleSpringModel.quat_from_axisangle(v,π/2)

bead = MicrotubuleSpringModel.Bead(
    BeadPos(0,0,0), 
    q,
    false, 
    false,
    (0,0), 0, 0
)

torque = MVector{3, Float64}(0,1,0)
q_τ = quat(0, torque...)
bead.q = sign(bead.q)
bead.q = bead.q - 0.5*q_τ * bead.q * dt
bead.q = sign(bead.q)

MicrotubuleSpringModel.to_axis_angle(bead.q)


# q = MicrotubuleSpringModel.quat_from_axisangle([0,1,0],π/2)

# v_ = MicrotubuleSpringModel.transform_orientation(v, q)

GLMakie.activate!()
GLMakie.closeall()
scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
vecs = normalize([imag_part(bead.q)...])
arrows!(scene, [Point3f(bead.x)], Vector{Vec{3, Float32}}([vecs]), linewidth=0.2, color=:white, arrowsize=0.3)
arrows!(scene, [Point3f(bead.x)], Vector{Vec{3, Float32}}([v]), linewidth=0.2, color=:red, arrowsize=0.3)
arrows!(scene, [Point3f(bead.x)], Vector{Vec{3, Float32}}([torque]), linewidth=0.2, color=:blue, arrowsize=0.3)
GLMakie.scale!(scene, 0.05, 0.05, 0.05)
mesh!(scene, Sphere(Point3f(bead.x), 0.2), color=MicrotubuleSpringModel.COLORS[(bead.α, bead.kinesin)], shininess=32.0)
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
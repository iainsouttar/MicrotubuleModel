

using GLMakie
using Parameters

GLMakie.activate!()

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)

a = 2.55
idx = 9*13
@unpack north, south, east, west = bead_info[idx]

bonds = [north, east, south, west]
colors = MicrotubuleSpringModel.NATURE.colors[4:7]


s = Scene(scene, camera=scene.camera)
plots = [mesh!(s, Sphere(Point3f(beads[b].x), a), color=c) for (b,c) in zip(bonds,colors) if b!=0]

scene


using GLMakie
GLMakie.activate!(ssao=true)
GLMakie.closeall() # close any open screen

fig = Figure()
ssao = Makie.SSAO(radius = 5.0, blur = 3)
ax = LScene(fig[1, 1], scenekw = (ssao=ssao,))
# SSAO attributes are per scene
ax.scene.ssao.bias[] = 0.025

box = Rect3(Point3f(-0.5), Vec3f(1))
positions = [Point3f(x, y, rand()) for x in -5:5 for y in -5:5]
meshscatter!(ax, positions, marker=box, markersize=1, color=:lightblue, ssao=true)
fig

radiance = 50000
lights = [
    PointLight(Vec3f(50, 0, 200), RGBf(radiance, radiance, radiance*1.1)),
]

fig = Figure()
ssao = Makie.SSAO(radius = 5.0, blur = 3)
ax = LScene(fig[1, 1], scenekw = (ssao=ssao,lights=lights))
scene.ssao.bias[] = 0.025
#mesh!(ax, Sphere(Point3f(0,0,0), a/2), color=:red, shininess=1.0, ssao=true)

colors = colorschemes[:nature].colors
COLORS = Dict(
    (true, false) => colors[1], 
    (false, false) => colors[2],
    (true, true) => colors[3],
    (false, true) => colors[4]
)
for bead in lattice
    mesh!(ax, Sphere(Point3f(bead.x), a/2), color=COLORS[(bead.α, bead.kinesin)], shininess=32.0)
end


fig



########################################################################

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info, dirs)

using MicrotubuleSpringModel: NATURE

c = [NATURE.colors[3],NATURE.colors[4],NATURE.colors[5],NATURE.colors[6]]


F = zeros(Float64, (3, lastindex(lattice)))
torque = similar(F)

MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([2*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)

scene

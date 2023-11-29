

using GLMakie

GLMakie.activate!()

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")

conf = set_bond_angles(conf)

lattice, dirs = MicrotubuleSpringModel.initialise(conf)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(lattice)

a = 2.55
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





GLMakie.activate!()
GLMakie.closeall()
scene = plot(lattice, dirs)

using MicrotubuleSpringModel: NATURE

c = [NATURE.colors[3],NATURE.colors[4],NATURE.colors[5],NATURE.colors[6]]

for (i,b) in enumerate(lattice)
    lat, long = b.lat_nn, b.long_nn
    intra = b.intra_nn
    bonds = [intra,lat[2],long, lat[1]]
    for (j,(bond,dir)) in enumerate(zip(bonds,eachcol(dirs[b.α])))
        if bond != 0
            arrows!(scene, [Point3f(b.x)], Vector{Vec{3, Float32}}([5.5*nhat[:,j, i]]), linewidth=0.35, arrowsize=0.5, color=c[j])
        end
    end
end


F = zeros(Float64, (3, lastindex(lattice)))
torque = similar(F)

MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, lattice, dirs, conf.spring_consts)

xs = [Point3f(b.x) for b in lattice]
vs = Vector{Vec{3, Float32}}([2*t for t in eachcol(torque)])
arrows!(scene, xs, vs, linewidth=0.3, arrowsize=0.5, color=:white)

scene

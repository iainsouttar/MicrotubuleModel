

function CairoMakie.Point3f(b::BeadPos)
    return Point3f(b.x, b.y, b.z)
end

function CairoMakie.Point2f(b::BeadPos)
    return Point2f(b.x, b.y)
end

scene = Scene(resolution=(1200,900))
cam3d!(scene)


s1 = Scene(scene, camera=scene.camera)
sphere_plot = mesh!(s1, Sphere(Point3f(b1), 0.5), color=:red)
s2 = Scene(scene, camera=scene.camera)
sphere_plot = mesh!(s1, Sphere(Point3f(b2), 0.5), color=:blue)
s3 = Scene(scene, camera=scene.camera)
sphere_plot = mesh!(s1, Sphere(Point3f(b3), 0.5), color=:yellow)


scale!(scene, 0.5, 0.5, 0.5)
#rotate!(scene, Vec3f(1, 0, 0), 0.3) # 0.5 rad around the y axis
scene


using GLMakie

GLMakie.activate!()

function plot_lattice(lattice)
    scene = Scene(resolution=(1200,900))
    cam3d!(scene)

    COLORS = Dict(true => :blue, false => :red)

    for bead in lattice
        s = Scene(scene, camera=scene.camera)
        sphere_plot = mesh!(s, Sphere(Point3f(bead.x), a/2), color=COLORS[bead.α])
    end


    scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end

plot_lattice(lattice)
scene


scene = Scene(resolution=(1200,900))
cam3d!(scene)

for pos in positions
    s = Scene(scene, camera=scene.camera)
    sphere_plot = mesh!(s, Sphere(Point3f(pos), a/2))
end


scale!(scene, 0.05, 0.05, 0.05)
center!(scene)
return scene

lattice = create_lattice(7, a, δx; S=3, N=13)
scene = plot_lattice(lattice)

idx = 13
lat, long = neighbours(idx, length(lattice))

s = Scene(scene, camera=scene.camera)
lat_plots = [mesh!(s, Sphere(Point3f(lattice[b].x), a/2), color=:black) for b in lat]
long_plots = [mesh!(s, Sphere(Point3f(lattice[b].x), a/2), color=:yellow) for b in long]
mesh!(s, Sphere(Point3f(lattice[idx].x), a/2), color=:pink)

scene

Ntot = length(lattice)

pos = [b.x for b in lattice]

distances = reshape([dot(x1-x2,x1-x2) for x1 in pos for x2 in pos],(Ntot,Ntot))
distances[distances .== 0] .= Inf

for i in 1:Ntot
    lat, long = neighbours(i, Ntot)
    neigh = Set((lat...,long...))
    nn = Set(sortperm(distances[:,i])[1:4])
    @show i, setdiff(neigh,nn)
end

neigh = Set((lat...,long...))
nn = Set(sortperm(distances[:,13])[1:4])

setdiffs = 0
for i in 1:Ntot
    lat, long = neighbours(i, Ntot)
    neigh = Set((lat...,long...))
    nn = Set(sortperm(distances[:,i])[1:4])
    setdiffs += length(setdiff(neigh,nn))
end
setdiffs == 0
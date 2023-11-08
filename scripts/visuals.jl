

using GLMakie

GLMakie.activate!()


lattice = create_lattice(7, a, δx; S=3, N=13)
scene = plot_lattice(lattice)

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

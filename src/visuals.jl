
function GLMakie.plot(lattice; a=4.05)
    scene = Scene(resolution=(1200,900))
    plot!(scene, lattice; a=a)
    return scene
end


function GLMakie.plot!(scene, lattice; a=4.05)
    colors = colorschemes[:seaborn_bright].colors
    COLORS = Dict(true => colors[5], false => colors[3])

    for bead in lattice
        s = Scene(scene, camera=scene.camera)
        sphere_plot = mesh!(s, Sphere(Point3f(bead.x), a/2), color=COLORS[bead.α])
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end
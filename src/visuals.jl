
colorschemes[:nature] = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

function GLMakie.plot(lattice; a=4.05)
    scene = Scene(resolution=(1200,900))
    plot!(scene, lattice; a=a)
    return scene
end


function GLMakie.plot!(scene, lattice; a=4.05)
    colors = colorschemes[:nature].colors
    COLORS = Dict(
        (true, false) => colors[1], 
        (false, false) => colors[2],
        (true, true) => colors[3],
        (false, true) => colors[4]
    )
    s = Scene(scene, camera=scene.camera)
    for bead in lattice
        mesh!(s, Sphere(Point3f(bead.x), a/2), color=COLORS[(bead.α, bead.kinesin)], shininess=32.0)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end
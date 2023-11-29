
ColorSchemes.colorschemes[:nature] = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

NATURE = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

COLORS = Dict(
    (true, false) => NATURE.colors[1], 
    (false, false) => NATURE.colors[2],
    (true, true) => NATURE.colors[3],
    (false, true) => NATURE.colors[4]
)

function GLMakie.plot(lattice::Vector{Bead}; a=4.05)
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    plot!(scene, lattice; a=a)
    return scene
end

function GLMakie.plot(lattice::Vector{Bead}, dirs; a=4.05, l=4.0)
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    plot!(scene, lattice, dirs; a=a, l=l)
    return scene
end


function GLMakie.plot!(scene, lattice::Vector{Bead}; a=4.05)
    for b in lattice
        mesh!(scene, Sphere(Point3f(b.x), a/4), color=COLORS[(b.α, b.kinesin)], shininess=32.0)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end

function GLMakie.plot!(scene, lattice::Vector{Bead}, dirs; a=4.05, l=4.0)
    c = [NATURE.colors[3],NATURE.colors[4],NATURE.colors[5],NATURE.colors[6]]

    for b in lattice
        for (i,bond) in enumerate(eachcol(dirs[b.α]))
            v = MicrotubuleSpringModel.transform_orientation(bond,b.q)
            v_ = l*normalize([imag_part(v)...])
            arrows!(scene, [Point3f(b.x)], Vector{Vec{3, Float32}}([v_]), linewidth=0.2, color=c[i], arrowsize=0.3)
        end
    end

    plot!(scene, lattice; a=a)
    return scene
end
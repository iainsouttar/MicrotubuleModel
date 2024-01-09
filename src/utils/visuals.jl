
ColorSchemes.colorschemes[:nature] = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

NATURE = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

COLORS = Dict(
    (true, false) => NATURE.colors[1], 
    (false, false) => NATURE.colors[2],
    (true, true) => NATURE.colors[3],
    (false, true) => NATURE.colors[4]
)

function GLMakie.plot(lattice::Vector{Bead}, info::Vector{BeadPars}; a=8.05)
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    plot!(scene, lattice, info; a=a)
    return scene
end

function GLMakie.plot(lattice::Vector{Bead}, info::Vector{BeadPars}, dirs; a=4.05, l=4.0)
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    plot!(scene, lattice, info, dirs; a=a, l=l)
    return scene
end


function GLMakie.plot!(scene, lattice::Vector{Bead}, info::Vector{BeadPars}; a=8.05)
    for (b,b_) in zip(lattice, info)
        mesh!(scene, Sphere(Point3f(b.x), a/4), color=COLORS[(b_.α, b.kinesin)], shininess=32.0)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end

function GLMakie.plot!(scene, lattice::Vector{Bead}, info::Vector{BeadPars}, dirs; a=4.05, l=4.0)
    c = [NATURE.colors[3],NATURE.colors[4],NATURE.colors[5],NATURE.colors[6]]

    for (b,b_) in zip(lattice, info)
        for (i,bond) in enumerate(eachcol(dirs[b_.α]))
            v = MicrotubuleSpringModel.transform_orientation(bond,b.q)
            v_ = l*normalize([imag_part(v)...])
            arrows!(scene, [Point3f(b.x)], Vector{Vec{3, Float32}}([v_]), linewidth=0.2, color=c[i], arrowsize=0.3)
        end
    end

    plot!(scene, lattice, info; a=a)
    return scene
end


function plot_individual!(scene, lattice::Vector{Bead}, info::Vector{BeadPars}; a=8.05)
    beads = []
    for (b,b_) in zip(lattice, info)
        bead = mesh!(scene, Sphere(Point3f(b.x), a/2), color=COLORS[(b_.α, b.kinesin)], shininess=32.0)
        push!(beads, bead)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return beads
end

function plot_flat!(ax::Axis, beads::Vector{Bead}, bead_info::Vector{BeadPars}; markersize=40)
    pts = [Point2f(b.x[1],b.x[3]) for b in beads]
    color = [COLORS[(b_.α, b.kinesin)] for (b,b_) in zip(beads, bead_info)]
    return scatter!(ax, pts, color=color, marker=:circle, markersize=markersize)
end



function plot_E!(ax, time, E)
    labels = [L"E_{lat}^r", L"E_{long}^r", L"E_{in}^r", L"E_{lat}^\theta", L"E_{long}^\theta", L"E_{in}^\theta"]
    E_band = cumsum(E,dims=1)
    E_tot = E_band[end,:]
    lines!(ax,time, E_tot, color=:black, linewidth=4)
    band!(ax, time, 0, E_band[1,:], color=NATURE.colors[1], label=labels[1])
    for i in 1:5
        lines!(ax,time, E_band[i,:], color=:black, linewidth=2)
        band!(ax, time, E_band[i,:], E_band[i+1,:], color=MicrotubuleSpringModel.NATURE.colors[i+1], label=labels[i+1])
    end
    xlims!(0,time[end])
    ylims!(low=0.0)
    ax.xlabel = "Iteration number"
    ax.ylabel = "Total Energy"
    axislegend(ax, position=:rt, nbanks = 2)
    return ax
end
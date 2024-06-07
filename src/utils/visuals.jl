
ColorSchemes.colorschemes[:nature] = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

COLORSCHEME = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

BEADCOLORS = Dict(
    (true, false) => COLORSCHEME.colors[1], 
    (false, false) => COLORSCHEME.colors[2],
    (true, true) => COLORSCHEME.colors[3],
    (false, true) => COLORSCHEME.colors[4]
)

function GLMakie.plot(lattice::Lattice, info::Vector{BeadPars}; a=4.05, l=4.0)
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    plot!(scene, lattice, info; a=a, l=l)
    return scene
end



function GLMakie.plot!(scene, lattice::Lattice, info::Vector{BeadPars}; a=4.05, l=4.0)
    c = [COLORSCHEME.colors[3],COLORSCHEME.colors[4],COLORSCHEME.colors[5],COLORSCHEME.colors[6]]

    for (x, q, b_) in zip(lattice.x, lattice.q, info)
        for (i,bond) in enumerate(b_.directions)
            v = MicrotubuleSpringModel.transform_orientation(bond,q)
            v_ = l*normalize([imag_part(v)...])
            arrows!(scene, [Point3f(x)], Vector{Vec{3, Float32}}([v_]), linewidth=0.2, color=c[i], arrowsize=0.3)
        end
    end

    for (x,kin, b_) in zip(lattice.x, lattice.kinesin, info)
        mesh!(scene, Sphere(Point3f(x), a/8), color=BEADCOLORS[(b_.α, kin)], shininess=32.0)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return scene
end


function plot_individual!(scene, lattice::Lattice, info::Vector{BeadPars}; a=8.05)
    beads = []
    for (x, kin, b_) in zip(lattice.x, lattice.kinesin, info)
        bead = mesh!(scene, Sphere(Point3f(x), a/4), color=BEADCOLORS[(b_.α, kin)], shininess=32.0)
        push!(beads, bead)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return beads
end


function plot_individual!(scene, x; a=8.05)
    beads = []
    colors = Dict(true=>COLORSCHEME.colors[1], false=>COLORSCHEME.colors[2])
    for (i,x_i) in enumerate(x)
        bead = mesh!(scene, Sphere(Point3f(x_i), a/2), color=colors[iseven(i÷13)], shininess=32.0)
        push!(beads, bead)
    end

    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    return beads
end


function plot_flat!(ax::Axis, lattice::Lattice, bead_info::Vector{BeadPars}; markersize=40)
    pts = Vector{Point2f}([Point2f(x[1],x[3]) for x in lattice.x])
    color = [BEADCOLORS[(b_.α, kin)] for (kin,b_) in zip(lattice.kinesin, bead_info)]
    return scatter!(ax, pts, color=color, marker=:circle, markersize=markersize)
end



function plot_energy!(ax, time, E)
    labels = [L"E_{lat}^r", L"E_{long}^r", L"E_{in}^r", L"E_{lat}^\theta", L"E_{long}^\theta", L"E_{in}^\theta", L"Torsion"]
    E_band = cumsum(E,dims=1)
    E_tot = E_band[end,:]
    lines!(ax,time, E_tot, color=:black, linewidth=4)
    band!(ax, time, 0, E_band[1,:], color=COLORSCHEME.colors[1], label=labels[1])
    for i in 1:6
        lines!(ax,time, E_band[i,:], color=:black, linewidth=2)
        band!(ax, time, E_band[i,:], E_band[i+1,:], color=MicrotubuleSpringModel.COLORSCHEME.colors[i+1], label=labels[i+1])
    end
    xlims!(0,time[end])
    ylims!(low=0.0)
    ax.xlabel = "Iteration number"
    ax.ylabel = "Total Energy"
    axislegend(ax, position=:rt, nbanks = 2)
    return ax
end
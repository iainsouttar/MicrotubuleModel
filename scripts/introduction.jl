#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


Nt = 1_000_000
stp = 100
time = collect(0:stp:Nt)



conf = from_toml(MicrotubuleConfig, "config/eulerMT142.toml")
#conf = set_bond_angles(conf)

lattice, bead_info = initialise(conf)

E = zeros((7,length(time)))
cosangles = zeros(length(time), length(lattice.x))
lngth = zeros(length(time))
E[:,1], cosangles[1, :], lngth[1] = total_energy(lattice, bead_info)

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        E[:,i÷stp+1], cosangles[i÷stp+1, :], lngth[i÷stp+1] = total_energy(lattice, bead_info)
    end
end

function plot_2d(lattice, bead_info)
    CairoMakie.activate!()
    f = Figure(resolution=(1000,1500))
    ax = Axis(f[1,1], aspect=DataAspect())
    plot_flat!(ax, lattice, bead_info, markersize=60)
    hidedecorations!(ax)
    hidespines!(ax)
    return f
end

function plot_3d(lattice, bead_info)
    GLMakie.activate!()
    GLMakie.closeall()
    scene = plot(lattice, bead_info)
    return scene
end
CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
plot_energy!(ax, time, E)
f

save("plots/energyMTprefnoT.png", f)


CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, time, lngth)
f
save("plots/MTlengthprefnoT.png", f)

f = plot_3d(lattice, bead_info)
f
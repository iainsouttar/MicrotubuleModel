"""
Find the equilibrium position for the protofilament and plot the energy over time as it equilibriates.

"""

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleConfig, "config/protofilamentv2.toml")
#conf = set_bond_angles(conf)

lattice, bead_info = initialise(conf)

Nt = 100000
stp = 1
time = collect(0:stp:Nt)
E = zeros((6,length(time)))
E[:,1] = total_energy(lattice, bead_info)

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        E[:,i÷stp+1] = total_energy(lattice, bead_info)
    end
end

#########################################################
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



scene = plot_3d(lattice, bead_info)
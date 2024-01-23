"""
Find the equilibrium position for the lattice and plot the energy over time as it equilibriates.

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

conf = from_toml(MicrotubuleConfig, "config/equilibrium.toml")
#conf = set_bond_angles(conf)

lattice, bead_info = MicrotubuleSpringModel.initialise(conf)

Nt = 200
step = 1
time = collect(0:step:Nt)
E = zeros((6,length(time)))
E[:,1] = total_energy(lattice, bead_info)

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % step == 0
        E[:,i÷step+1] = total_energy(lattice, bead_info)
    end
end

###########################################################

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
plot_energy!(ax, time, E)
f

scene = plot_3d(lattice, bead_info)
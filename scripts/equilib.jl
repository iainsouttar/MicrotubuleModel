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
conf = set_bond_angles(conf)

beads, bead_info = MicrotubuleSpringModel.initialise(conf)

beads_cpy = deepcopy(beads)
# for (i,b) in enumerate(beads)
#     if i % 13 ∈ (4,5,6,7,8,9)
#         b.kinesin = true
#     end
# end

Nt = 10
step = 10
time = collect(0:step:Nt)
E = zeros(length(time))
E[1] = total_energy(beads, bead_info, conf.spring_consts.K)

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, conf, conf.iter_pars)
    if i % step == 0
        E[i÷step+1] = total_energy(beads, bead_info, conf.spring_consts.K)
    end
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, time, E)
f

beads, bead_info, dirs, E, time = burnin(conf, 1e4, show_E=true)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)
scene

###########################################################

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
MicrotubuleSpringModel.plot_E!(ax, time, E)
f
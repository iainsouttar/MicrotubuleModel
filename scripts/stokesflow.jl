#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stokes_flow.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

###################################################################

Nt = 50_000

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
end

CairoMakie.activate!()
f = Figure(resolution=(1000,1500))
ax = Axis(f[1,1], aspect=DataAspect())
plot!(ax, beads, bead_info)
hidedecorations!(ax)
hidespines!(ax)
f
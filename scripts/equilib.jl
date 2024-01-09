#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")
conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

# for (i,b) in enumerate(beads)
#     if i % 13 ∈ (4,5,6,7,8,9)
#         b.kinesin = true
#     end
# end

Nt = Int(1e6)

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
#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using ProgressMeter

FILENAME = "anim3d"

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stochastic.toml")

conf = set_bond_angles(conf)

lattice, bead_info = MicrotubuleSpringModel.initialise(conf)

Nt = 500
stp = 10
time = 0:stp:Nt

GLMakie.activate!()
GLMakie.closeall()
scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
bead_plots = plot_individual!(scene, beads, bead_info)
scene

p = Progress(Nt)
Makie.record(scene, "figures/$FILENAME.mp4", 1:Nt) do i
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    for (n, b) in enumerate(beads)
        translate!(bead_plots[n], b.x...)
    end
    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    next!(p)
end

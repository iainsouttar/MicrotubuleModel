#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using ProgressMeter

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stochastic.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

# for b in beads
#     b.x += 0.5*randn(3)
# end

Nt = 500
step = 10
time = 0:step:Nt

E = zeros((6, length(time)))

E[:,1] = total_energy(beads, bead_info, dirs, conf.spring_consts)

GLMakie.activate!()
GLMakie.closeall()
scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
bead_plots = plot_individual!(scene, beads, bead_info)
scene

p = Progress(Nt)
Makie.record(scene, "figures/anim2.mp4", 1:Nt) do i
    MicrotubuleSpringModel.iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    #@info deltax[2], deltax[1]
    for (n, b) in enumerate(beads)
        translate!(bead_plots[n], b.x...)
    end
    GLMakie.scale!(scene, 0.05, 0.05, 0.05)
    center!(scene)
    E[:,i÷step+1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
    next!(p)
end



GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)
scene

GLMakie.activate!()
GLMakie.closeall()
scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
cam3d!(scene)
bead1 = mesh!(scene, Sphere(Point3f(0,0,0), 2.0), color=MicrotubuleSpringModel.NATURE.colors[1], shininess=32.0)
bead2 = mesh!(scene, Sphere(Point3f(4,0,0), 2.0), color=MicrotubuleSpringModel.NATURE.colors[2], shininess=32.0)
GLMakie.scale!(scene, 0.05, 0.05, 0.05)
scene

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using CairoMakie

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/youngs_modulus.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

Nt = 2000
step = 10
time = 0:step:Nt

beads_obs = Observable(beads)
pts = @lift([Point2f(b.x[1],b.x[3]) for b in $(beads_obs)])
color = [MicrotubuleSpringModel.COLORS[(b_.α, b.kinesin)] for (b,b_) in zip(beads, bead_info)]

f = Figure(resolution=(1500,1500), backgroundcolor=colorant"#111111")
ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
#limits!(-14,50,-5,100)
limits!(-15,15,-5,55)
hidedecorations!(ax)
hidespines!(ax)
scatter!(ax, pts, color=color, marker=:circle, markersize=40)
f

p = Progress(Nt÷step+1)
frames = range(1, Nt÷step+1, step=1)
record(f, "figures/anim_flat_YM.mp4", frames, framerate=25) do i
    for t in 1:step
        MicrotubuleSpringModel.iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    end

    beads_obs[] = beads
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

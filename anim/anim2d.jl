#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using CairoMakie
CairoMakie.activate!()

FILENAME = "anim2d"

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/youngs_modulus.toml")
conf = set_bond_angles(conf)
beads, bead_info = MicrotubuleSpringModel.initialise(conf)

Nt = 2000
stp = 10
time = 0:stp:Nt

beads_obs = Observable(beads)
pts = @lift([Point2f(b.x[1],b.x[3]) for b in $(beads_obs)])

# 4 colours dependent on alpha/beta and kinesin bound or not
color = [MicrotubuleSpringModel.COLORS[(b_.α, b.kinesin)] for (b,b_) in zip(beads, bead_info)]

# plot original configuration
f = Figure(resolution=(1500,1500), backgroundcolor=colorant"#111111")
ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
limits!(-15,15,-5,55)
hidedecorations!(ax)
hidespines!(ax)
scatter!(ax, pts, color=color, marker=:circle, markersize=40)
f

# simulate Nt iterations and visualise every "stp" iterations 
p = Progress(Nt÷stp+1)
frames = range(1, Nt÷stp+1, step=1)
record(f, "figures/$FILENAME.mp4", frames, framerate=25) do i
    for t in 1:stp
        iterate!(beads, bead_info, conf, conf.iter_pars)
    end

    beads_obs[] = beads
    next!(p)
end

using CairoMakie
CairoMakie.activate!()

conf = from_toml(MicrotubuleConfig, "config/protofilament.toml")
conf = set_bond_angles(conf)
beads, bead_info = MicrotubuleSpringModel.initialise(conf)

Nt = 1_000
stp = 200
time = 0:stp:Nt

for i in 1:Nt
    iterate!(beads, bead_info, conf_burnin, conf_burnin.iter_pars)
end

beads_obs = Observable(beads)
pts = @lift(Vector{Point2f}([Point2f(x[1],x[3]) for x in $(beads_obs).x]))

# 4 colours dependent on alpha/beta and kinesin bound or not
color = [MicrotubuleSpringModel.COLORS[(b_.α, kin)] for (kin,b_) in zip(beads.kinesin, bead_info)]

# plot original configuration
f = Figure(resolution=(1500,1500), backgroundcolor=colorant"#111111")
ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
#limits!(-15,25,-5,40)
limits!(-15,15,-10,10*4.05+10)
hidedecorations!(ax)
hidespines!(ax)
scatter!(ax, pts, color=color, marker=:circle, markersize=100)
f



# simulate Nt iterations and visualise every "stp" iterations 
p = Progress(Nt÷stp+1)
frames = range(1, Nt÷stp+1, step=1)
record(f, "figures/protofilament-det.mp4", frames, framerate=25) do i
    for t in 1:stp
        iterate!(beads, bead_info, conf, conf.iter_pars)
    end

    beads_obs[] = beads
    next!(p)
end
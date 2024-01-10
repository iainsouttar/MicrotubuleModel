if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


function main!(F, conf, num_rings, step, Nt)
    conf = @set conf.lattice.num_rings = num_rings

    beads, bead_info = burnin(conf, 500)

    L0 = microtubule_length(beads, conf.lattice)
    N = conf.lattice.N
    Ntot = lastindex(beads)
    time = collect(0:step:Nt)
    deflection = zeros(length(time))

    original = deepcopy(beads[Ntot-N:end])
    
    #conf = @set conf.external_force = ym_conf

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf, conf.iter_pars)
        if i % step == 0
            deflection[i÷step+1] = deflection_end(beads, original)
        end
    end
    return stiffness(F, L0, deflection[end]), L0
end

Nt = 100_000
step = 200
time = 0:step:Nt

num_rings = collect(10:10:50)
L0 = zeros(length(num_rings))
stiff = zeros(length(num_rings))

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/bending_stiffness.toml")
conf = set_bond_angles(conf)

for (i, n) in enumerate(num_rings)
    stiff[i], L0[i] = main!(0.05, conf, n, step, Nt)
end


CairoMakie.activate!()
f = Figure(resolution=(1000,1500))
ax = Axis(f[1,1], aspect=DataAspect())
plot_flat!(ax, beads, bead_info)
hidedecorations!(ax)
hidespines!(ax)
f


f = Figure(resolution=(1000,700))
ax = Axis(f[1,1])
scatter!(ax, L0, stiff)
f
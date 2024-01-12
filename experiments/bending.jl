using Logging
using Parameters: @unpack
using ProgressMeter
using Configurations
using DataFrames


#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


function main!(F, conf, num_rings, step, Nt)
    conf = @set conf.lattice.num_rings = num_rings
    conf_burnin = deepcopy(conf)
    conf_burnin = @set conf_burnin.external_force = MicrotubuleSpringModel.NoExternalForce()
    beads, bead_info = burnin(conf_burnin, 500)

    L0 = microtubule_length(beads, conf.lattice)
    N = conf.lattice.N
    Ntot = length(beads)
    time = collect(0:step:Nt)
    deflection = zeros(length(time))

    original = deepcopy(beads.x[Ntot-N:end])

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf, conf.iter_pars)
        if i % step == 0
            deflection[i÷step+1] = deflection_end(beads.x, original)
        end
    end
    return stiffness(13*F, L0, deflection[end]), L0
end

Nt = 50_000
step = 2_000
time = 0:step:Nt

num_rings = collect(10:10:50)
L0 = zeros(length(num_rings))
stiff = zeros(length(num_rings))

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/bending_stiffness.toml")
conf = set_bond_angles(conf)

for (i, n) in enumerate(num_rings)
    stiff[i], L0[i] = main!(0.01, conf, n, step, Nt)
    @info stiff[i], L0[i]
end

df = DataFrame(length=L0, stiffness=stiff)
save_to_csv("bending_stiffness.csv", df)

# CairoMakie.activate!()
# f = Figure(resolution=(1000,1500))
# ax = Axis(f[1,1], aspect=DataAspect())
# plot_flat!(ax, beads, bead_info)
# hidedecorations!(ax)
# hidespines!(ax)
# f


# f = Figure(resolution=(1000,700))
# ax = Axis(f[1,1])
# scatter!(ax, L0, stiff)
# f
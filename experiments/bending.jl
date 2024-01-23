"""
Test the bending stiffness of microtubules of various lengths.

Finds the equilibrium length and position of the MT.
Applies a force to the free end perpendicular to MT.
Measures displacement of the end ring once new equilib has been found.

"""

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


function main!(F, conf, num_rings, Nt)
    conf = @set conf.lattice.num_rings = num_rings
    
    # first find equilib position of lattice
    conf_burnin = deepcopy(conf)
    conf_burnin = @set conf_burnin.external_force = MicrotubuleSpringModel.NoExternalForce()
    beads, bead_info = burnin(conf_burnin, 500)

    L0 = microtubule_length(beads, conf.lattice)
    N = conf.lattice.N
    Ntot = length(beads)
    original = deepcopy(beads.x[Ntot-N:end])

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf, conf.iter_pars)
    end

    return stiffness(13*F, L0, deflection_end(beads.x, original)), L0
end

Nt = 1_000
time = 0:stp:Nt
filename = "bending_stiffness_test.csv"

df = DataFrame(length=[0], stiffness=[0])
save_to_csv(filename, df)

num_rings = collect(10:10:50)
L0 = zeros(length(num_rings))
stiff = zeros(length(num_rings))

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/bending_stiffness.toml")
conf = set_bond_angles(conf)

for (i, n) in enumerate(num_rings)
    stiff[i], L0[i] = main!(0.01, conf, n, Nt)
    @info stiff[i], L0[i]
    save_to_csv(filename, df, append=true, header=false)
end


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

using DataFrames

data, header = readdlm("results/processed/bending_stiffness2.csv", ',', header=true)

df = DataFrame(data, vec(header))


CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1],
        #ylabel="Strain",
        #xlabel="Stress (GPa)"
)
scatter!(ax, df.length, df.stiffness)
#limits!(ax,0.0, stress[end]*1.02 ,0.0, strain[end]*1.02)
f
#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

function main!(beads, bead_info, F, conf, dirs, step, Nt, L0)
    ext = zeros(Nt÷step+1)
    conf = @set conf.external_force = MicrotubuleSpringModel.YoungsModulusTest(F, 13)
    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
        if i % step == 0
            ext[i÷step+1] = microtubule_length(beads, conf.lattice) - L0
        end
    end
    return ext
end

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/youngs_modulus.toml")
conf = set_bond_angles(conf)
beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

microtubule_length(beads, conf.lattice) ≈ (conf.lattice.num_rings-1) *conf.lattice.a

ym_conf = deepcopy(conf.external_force)
conf = @set conf.external_force = MicrotubuleSpringModel.NoExternalForce()

Nt = 500
step = 1
time = 0:step:Nt
E = zeros((6,Nt÷step+1))

E[:,1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    if i % step == 0
        E[:,i÷step+1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
    end
end

beads_cpy = deepcopy(beads)
bead_info_cpy = deepcopy(bead_info)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
MicrotubuleSpringModel.plot_E!(ax, time, E)
f

###################################################################

conf = @set conf.external_force = ym_conf

L0 = microtubule_length(beads, conf.lattice)

forces = 0.1:0.2:0.6
Nt = 10_000
step = 50
time = 0:step:Nt

stress = forces ./ SA(10.61,2)
strain = zeros(length(forces))

for (i,F) in enumerate(forces)
    ext = main!(beads, bead_info, F, conf, dirs, step, Nt, L0)
    strain[i] = ext[end] / L0
end

using DataFrames, GLM, StatsBase

data = DataFrame(X=stress, Y=strain)
sol = lm(@formula(Y ~ X), data)

YM = 1/ coef(sol)[2]
print(YM, " +/- " , stderror(sol)[2]*YM^2)

x = 0.0:0.00001:maximum(stress)*1.05
est = predict(sol, DataFrame(X=x), interval=:confidence)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1],
        ylabel="Strain",
        xlabel="Stress (GPa)"
)
band!(ax, x, est.lower, est.upper, color=MicrotubuleSpringModel.NATURE.colors[2])
lines!(ax, x, est.prediction, linewidth=3, color=:black)
scatter!(ax, data.X, data.Y)
limits!(ax,0.0, stress[end]*1.02 ,0.0, strain[end]*1.02)
f
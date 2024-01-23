"""
Test the Young's modulus of microtubules of various lengths.

Finds the unstretched equilibrium length of the MT.
Applies a force to the free end parallel to MT axis.
Measures change in length once new equilib has been found.

"""

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

function main!(beads, bead_info, F, conf, stp, Nt, L0)
    ext = zeros(Nt÷stp+1)
    conf = Setfield.@set conf.external_force = MicrotubuleSpringModel.YoungsModulusTest(F, 13)
    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, conf, conf.iter_pars)
        if i % stp == 0
            ext[i÷stp+1] = microtubule_length(beads, conf.lattice) - L0
        end
    end
    return ext
end

#####################################################

conf = from_toml(MicrotubuleConfig, "config/youngs_modulus.toml")
conf = set_bond_angles(conf)

beads, bead_info = burnin(conf, 2_000)

###################################################################

L0 = microtubule_length(beads, conf.lattice)

forces = 0.0:0.1:0.5
Nt = 50_000
stp = 50

stress = 13*forces ./ surface_area(10.61,2)
strain = zeros(length(forces))

for (i,F) in enumerate(forces)
    ext = main!(beads, bead_info, F, conf, stp, Nt, L0)
    @info F, ext[end]
    strain[i] = ext[end] / L0
end

####################################################################

using DataFrames, GLM, StatsBase

function fit_YM(data)
    sol = lm(@formula(Y ~ X), data)
    YM = 1/ coef(sol)[2]
    return sol, YM, stderror(sol)[2]*YM^2
end

data = DataFrame(X=stress, Y=strain)
sol, YM, err = fit_YM(data)

print(YM, " +/- " , err)

x = 0.0:0.00001:maximum(stress)*1.05
est = predict(sol, DataFrame(X=x), interval=:confidence)


CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1],
        ylabel="Strain",
        xlabel="Stress (GPa)"
)
band!(ax, x, est.lower, est.upper, color=COLORSCHEME.colors[2])
lines!(ax, x, est.prediction, linewidth=3, color=:black)
scatter!(ax, stress, strain)
limits!(ax,0.0, stress[end]*1.02 ,0.0, strain[end]*1.02)
f

save_to_csv("young-modulus-fit.csv", collect(x), [est.prediction, est.lower, est.upper])
save_to_csv("young-modulus-data.csv", stress, [strain])

using CSV 

est[!, "x"] = collect(x)
select!(est, :x, Not([:x]))
save_to_csv("young-modulus-fit-2.csv", est)

save_to_csv("young-modulus-data-2.csv", DataFrame(stress=stress, strain=strain))


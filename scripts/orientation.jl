#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/rotation.toml")
conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

# for (i,b) in enumerate(beads)
#     if i % 13 ∈ (4,5,6,7,8,9)
#         b.kinesin = true
#     end
# end

Nt = Int(1e6)
step = 1
time = 0:step:Nt
E = zeros((6, length(time)))

E[:,1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    if i % step == 0
        E[:,i÷step+1] = total_energy(beads, bead_info, dirs, conf.spring_consts)
    end
end

F = zeros(Float64, (3, lastindex(beads)))
torque = similar(F)
MicrotubuleSpringModel.eval_forces_and_torques!(F, torque, beads, bead_info, dirs, conf.spring_consts)

GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)
scene

############################################################

labels = [L"E_{lat}^r", L"E_{long}^r", L"E_{in}^r", L"E_{lat}^\theta", L"E_{long}^\theta", L"E_{in}^\theta"]

E_tot = vec(sum(E,dims=1))
idx = argmin(E_tot)
CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1], xlabel="Iteration number", ylabel="Total Energy")
lines!(ax,time, E_tot,color=:black, linewidth=4, label=L"E_T")
series!(ax,time,E, color=[MicrotubuleSpringModel.NATURE.colors...,:red],linewidth=4, labels=labels)
vlines!(ax, time[idx], 0.0,25)
xlims!(0,Nt/2)
ylims!(low=0.0)
axislegend(ax, position=:rt)
f

save("non-increasing-energy.png", f)

###########################################################

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
MicrotubuleSpringModel.plot_E!(ax, time, E)
f


k_ab = 3.5
k_bb = 1.0/0.4
k_a = 3.6
k_b = 3.8
K_ab = 1/(1/(2*k_a)+1/(2*k_b)+1/k_ab)

K_bb = 1/(1/k_b+1/k_bb)

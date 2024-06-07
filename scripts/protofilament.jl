"""
Find the equilibrium position for the protofilament and plot the energy over time as it equilibriates.

"""

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleConfig, "config/protofilamentv2.toml")
#conf = set_bond_angles(conf)
path = "results/raw"
filename = "fluctuations-free-test_startNotorsion.csv"

lattice, bead_info = initialise(conf)
Nt = 90000
stp = 100
time = collect(0:stp:Nt)
E = zeros((7,length(time)))
cosangles = zeros(length(time), length(lattice.x))
lngth = zeros(length(time))
E[:,1], cosangles[1, :], lngth[1] = total_energy(lattice, bead_info)

data = Matrix{Float64}(zeros(Float64, (length(lattice.x)*3,1)))
for i in 1:length(lattice.x)
    data[3*(i-1)+1:3*i,1] .= lattice.x[i]
end
open(path*"/"*filename, "w") do io
    writedlm(io, data', ',')
end

@showprogress for i in 1:Int(Nt/2)
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        E[:,i÷stp+1], cosangles[i÷stp+1, :],lngth[i÷stp+1] = total_energy(lattice, bead_info)
        MicrotubuleSpringModel.append_to_csv(filename, lattice.x)

    end
end

# l0_kin = conf.spring_consts.l0_in_kin
# k_kin = conf.spring_consts.k_in_kin
# K_in_kin = conf.spring_consts.K_in_kin
# for i in 1:length(lattice.x)
#     #print(i)
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i))
# end



@showprogress for i in Int(Nt/2 + 1):(Nt)
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        E[:,i÷stp+1], cosangles[i÷stp+1, :], lngth[i÷stp+1] = total_energy(lattice, bead_info)
        MicrotubuleSpringModel.append_to_csv(filename, lattice.x)
        #print(E[7,i÷stp+1])

    end
end


#########################################################
function plot_2d(lattice, bead_info)
    CairoMakie.activate!()
    f = Figure(resolution=(1000,1500))
    ax = Axis(f[1,1], aspect=DataAspect())
    plot_flat!(ax, lattice, bead_info, markersize=60)
    hidedecorations!(ax)
    hidespines!(ax)
    return f
end

function plot_3d(lattice, bead_info)
    GLMakie.activate!()
    GLMakie.closeall()
    scene = plot(lattice, bead_info)
    return scene
end


CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
plot_energy!(ax, time, E)
f
save("plots/zigzagtest1.png", f)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = f[1,1] = Axis(f)
for i in 1:length(lattice.x)
    lines!(ax, time, cosangles[:, i], label="$i")
end
f[1,2] = Legend(f, ax, framevisible = false)
f

#save("plots/zigzag_angles_19.png", f)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, time, lngth)
f
save("plots/zigzag_lengths_test1.png", f)


cutoff = 8
CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1], xlabel = "Monomer number", ylabel = "x value of position vector", title = "Zigzag")
lines!(ax,cutoff:(length(lattice.x)-cutoff+1) , getindex.(lattice.x,1)[cutoff:(length(lattice.x)-cutoff+1)])
lines!(ax, cutoff:(length(lattice.x)-cutoff+1), 4.05*sin(0.1)*ones(length(lattice.x)-2*cutoff+2))

f
save("plots/zigzag_x2.png", f)


zigzagKin = plot_3d(lattice, bead_info)
#save("plots/test.png", scene)

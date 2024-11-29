if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


Nt = 36000
stp = 10
time = collect(0:stp:Nt)
conf = from_toml(MicrotubuleConfig, "config/eulerMTkinesinLikestoch.toml")
reps = 1

function angle(a, b)
    return acosd(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
end


lattice, bead_info = initialise(conf)

# for i in 8:40
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+1), true)
# end

# for i in 8:40
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+3), true)
# end

# for i in 8:40
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+5), true)
# end

# for i in 8:40
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+7), true)
# end

# for i in 8:40
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+9), true)
# end

# for i in 25:30
#     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+11), true)
# end
strt = -3.0
nd = -3.0
num = 1
a1 = 10 .^collect(range(strt,nd,num))
a2 = 10 .^collect(range(strt,nd,num))

E = zeros((num,num,7,length(time)))
cosangles = zeros(length(time), length(lattice.x))
lngth = zeros(num,num,length(time))
energy_comps = zeros((num,num,length(time), 7, length(lattice.x)))
angles25 = zeros(num,num,(length(time)-1)*5, 12)
clusters = falses(num,num,length(time),length(lattice.x))
#lattice, bead_info = initialise(conf)
kin_index = 11*13+5
for (kout,trans_rate_out_kin) in enumerate([0.0])
    for (kin,trans_rate_in_kin) in enumerate([0.0])

        #lattice, bead_info = initialise(conf)


        E[kout,kin,:,1], cosangles[1, :], lngth[kout,kin,1], energy_comps[kout,kin,1,:,:] = total_energy(lattice, bead_info, conf.lattice.N,  conf.lattice.S)

        @showprogress for i in 1:Nt
            iterate!(lattice, bead_info, conf, conf.iter_pars,0, trans_rate_out_kin, trans_rate_in_kin, 0.0005) #1e-8

            # if i < 12000
            #     iterate!(lattice, bead_info, conf, conf.iter_pars,0, 0.0, 0.0, 0.0005) #1e-8
            # else
            #     iterate!(lattice, bead_info, conf, conf.iter_pars,0, trans_rate_out_kin, trans_rate_in_kin, 0.0005) #1e-8
            # end
                if i == 12000
                    
                for i in 0:47
                    attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+1), true)
                end

                # for i in 8:40
                #     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+3), true)
                # end

                # for i in 8:40
                #     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+5), true)
                # end

                for i in 0:47
                    attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(i*13+7), true)
                end
            end
            
            if i % stp == 0
                # if i==12000
                #     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(kin_index), true)
                #     print(bead_info[kin_index].directions, "\n", bead_info[kin_index+5].directions)
                # end
                # if i==24000
                #     remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(kin_index), true, false)
                # end
                    
                
                clusters[kout,kin,(i÷stp)+1,:] = lattice.kinesin[:,2]
                E[kout,kin,:,i÷stp+1], cosangles[i÷stp+1, :], lngth[kout,kin,i÷stp+1], energy_comps[kout,kin,i÷stp+1,:,:] = total_energy(lattice, bead_info, conf.lattice.N,  conf.lattice.S)
                for r in 1:5
                    ring = lattice.x[((25+r)*13+1):((25+r)*13+13)]
                    central2d = mean(ring)[1:2]

                    for j in 1:12
                        vec1 = ring[j][1:2] - central2d
                        vec2 = ring[(j+1)][1:2] - central2d
                        angles25[kout,kin,(i÷stp)+(length(time)-1)*(r-1),j] = angle(vec1,vec2)
                    end
                end
            end
        end
    end
end

CairoMakie.activate!()


fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1])
#ax = Axis(fig[1,1], aspect=AxisAspect(1))
hist!(ax, vcat([angles25[1,1,((length(time)-1)*r-1000):((length(time)-1)*r), [1]][:] for r in 1:5]...), bins = 100, normalization = :pdf, label = "right seam")
hist!(ax, vcat([angles25[1,1,((length(time)-1)*r-1000):((length(time)-1)*r), [1]][:] for r in 1:5]...), bins = 100, normalization = :pdf, label = "right seam")
hist!(ax, vcat([angles25[1,1,((length(time)-1)*r-1000):((length(time)-1)*r), 6:6][:] for r in 1:5]...), bins = 100, normalization = :pdf, label = "Nonseam")


ax.xlabel = "Wall angle between neighbouring tubulin"
ax.ylabel = "Frequency"
ax.title = "Wall angle with kinesin"
axislegend(ax, position =:rt)
#lines!(ax, angles[100,:])
display(fig)
#save("wall_angle.png", fig)

CairoMakie.activate!()


fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1])
#ax = Axis(fig[1,1], aspect=AxisAspect(1))
hist!(ax, (angles25[1,1,((length(time)-1)*1-1000):((length(time)-1)*1), 5]), bins = 100, normalization = :pdf, label = "right seam")
hist!(ax, (angles25[1,1,((length(time)-1)*1-1000):((length(time)-1)*1), 6]), bins = 100, normalization = :pdf, label = "right seam")
hist!(ax, (angles25[1,1,((length(time)-1)*1-1000):((length(time)-1)*1), 4]), bins = 100, normalization = :pdf, label = "right seam")
hist!(ax, (angles25[1,1,((length(time)-1)*1-1000):((length(time)-1)*1), 3]), bins = 100, normalization = :pdf, label = "right seam")

#ist!(ax, vcat([angles25[1,1,((length(time)-1)*r-1000):((length(time)-1)*r), [1]][:] for r in 1:5]...), bins = 100, normalization = :pdf, label = "right seam")
#hist!(ax, vcat([angles25[1,1,((length(time)-1)*r-1000):((length(time)-1)*r), 6:6][:] for r in 1:5]...), bins = 100, normalization = :pdf, label = "Nonseam")


ax.xlabel = "Wall angle between neighbouring tubulin"
ax.ylabel = "Frequency"
ax.title = "Wall angle with kinesin"
axislegend(ax, position =:rt)
#lines!(ax, angles[100,:])
display(fig)
#save("wall_angle.png", fig)


# fig = Figure(resolution = (1500, 1000))
# ax = Axis(fig[1,1], aspect=AxisAspect(1))

# scatter!(ax, ring2d)
# scatter!(ax, central2d[1],central2d[2])
# display(fig)


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

f1 = plot_3d(lattice, bead_info)
f1

CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])

ax.xlabel = "Time stamp"
ax.ylabel = "number of kinesin-bound-like tubulin"
ax.title = "cluster growth"

lines!(ax, vec(sum(clusters[1,1,:,:], dims = 2)))

#axislegend(ax, position=:lc)
display(fig)
#save("plots/cluster_growth_energy_kin14.png", fig)

CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
plot!(ax, clusters[1,1,:,:]) #5,9
display(fig)

CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))
energy_term_list = ["Lat len", "Lon len", "Int len", "Lat bend", "Lon bend", "Int bend"]

for i in 1:6
    lines!(ax, vec(E[1,1,i,:]), label = energy_term_list[i])
end
axislegend(ax, position=:rc)

display(fig)


fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))
proto_indices = zeros(Int64, (conf.lattice.num_rings, conf.lattice.N))
for i in 1:conf.lattice.N
    for j in 1:conf.lattice.num_rings
        proto_indices[j,i] = Int64((j-1)*conf.lattice.N + i)
    end
end

seam_right = proto_indices[:,conf.lattice.N]
seam_left = proto_indices[:,1]
rest = trues(length(lattice.x))
rest[seam_right] .=false
rest[seam_left] .=false
buffer = 4

rest[1:(buffer*conf.lattice.N-buffer-1)] .= false

rest[(conf.lattice.num_rings-buffer)*conf.lattice.N:length(lattice.x)]
times = 200:700

CairoMakie.activate!()

hist!(ax, vec(sum(energy_comps[1,1,times,:,:], dims = 2)[:,:,seam_right[buffer:(conf.lattice.num_rings-buffer)]]),normalization = :pdf, label = "right seam")
hist!(ax, vec(sum(energy_comps[1,1,times,:,:], dims = 2)[:,:,seam_left[buffer:(conf.lattice.num_rings-buffer)]]),normalization = :pdf, label = "left seam")
hist!(ax, vec(sum(energy_comps[1,1,times,:,:], dims = 2)[:,:,rest]), normalization = :pdf, label = "Non-seam")

display(fig)

CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))

energy_term_list = ["Lat len", "Lon len", "Int len", "Lat bend", "Lon bend", "Int bend"]
for bnd in (14*20+1):(14*20+14) #bead_info[kin_index].bonds
    for (p,energy_term) in enumerate(energy_term_list)
        #print(bnd)
        lines!(ax, vec(sum(energy_comps[1,1,:,p,bnd], dims = 2)), label = "$bnd,$energy_term")
    end
end
# lines!(ax, vec(sum(energy_comps[1,1,:,:,101], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,87], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,99], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,113], dims = 2)))

axislegend(ax, position=:rc)
display(fig)


CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))

for bnd in (14*20+1):(14*20+14) #bead_info[kin_index].bonds #bead_info[kin_index].bonds
        #print(bnd)
    lines!(ax, log.(vec(sum(energy_comps[1,1,:,:,bnd], dims = 2))), label = "$bnd")
end
# lines!(ax, vec(sum(energy_comps[1,1,:,:,101], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,87], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,99], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,113], dims = 2)))
#lines!(ax, log.(vec(sum(energy_comps[1,1,:,:,kin_index], dims = 2))), label = "$kin_index")


ax.ylabel = "Energy"
axislegend(ax, position=:rc)
display(fig)

#save("energy_time.png", fig)

CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))

for i in 1:13 #bead_info[kin_index].bonds #bead_info[kin_index].bonds
        #print(bnd)
    proto_energy = sum([energy_comps[1,1,:,:,j] for j in (i+8*13):13:(i+40*13)])
    print(size(proto_energy))
    lines!(ax, log.(vec(sum(proto_energy, dims = 2))), label = "$i")
end
# lines!(ax, vec(sum(energy_comps[1,1,:,:,101], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,87], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,99], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,113], dims = 2)))
#lines!(ax, log.(vec(sum(energy_comps[1,1,:,:,kin_index], dims = 2))), label = "$kin_index")


ax.ylabel = "Avg Energy across protofilament"
axislegend(ax, position=:rc)
display(fig)



CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))

for bnd in (kin_index .+ (-4:4)) #bead_info[kin_index].bonds
        #print(bnd)
    scatter!(ax,bnd-kin_index,  mean(vec(sum(energy_comps[1,1,650:1150,:,bnd], dims = 2))), label = "$bnd")
end
# lines!(ax, vec(sum(energy_comps[1,1,:,:,101], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,87], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,99], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,113], dims = 2)))
#lines!(ax, (vec(sum(energy_comps[1,1,:,:,kin_index], dims = 2))), label = "$kin_index")



#axislegend(ax, position=:rc)
display(fig)

#save("strain_along_helix_lateral.png", fig)

long_strain = [mean(vec(sum(energy_comps[1,1,650:1150,:,bnd], dims = 2))) for bnd in (kin_index .+ 13 .*(-6:6))]
lat_strain = [mean(vec(sum(energy_comps[1,1,650:1150,:,bnd], dims = 2))) for bnd in (kin_index .+ (-6:6))]
CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))

lines!(ax,-6:6, long_strain, label = "Longitudinal strain decay")
lines!(ax, -6:6, lat_strain, label = "Lateral strain decay")

# lines!(ax, vec(sum(energy_comps[1,1,:,:,101], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,87], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,99], dims = 2)))
# lines!(ax, vec(sum(energy_comps[1,1,:,:,113], dims = 2)))
#lines!(ax, (vec(sum(energy_comps[1,1,:,:,kin_index], dims = 2))), label = "$kin_index")
ax.xlabel = "Distance away from state 2 tubulin"
ax.ylabel = "Average energy in relaxed state"


axislegend(ax, position=:rc)
display(fig)

#save("strain_decay133.png", fig)

#counting seam

function rolling_average(y, window)
    return [mean(y[max(1, i - window + 1):i]) for i in 1:length(y)]
end
window_size = 100
CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1), xlabel = "timestep", ylabel = "Proportion of protofilament in kinesin-like")

#lines!(ax, rolling_average(vec(sum(clusters[1,1,:,seam_left[buffer:(conf.lattice.num_rings-buffer)]], dims = 2)./(conf.lattice.num_rings-2*buffer)), window_size), label = "1")
#lines!(ax, rolling_average(vec(sum(clusters[1,1,:,seam_right[buffer:(conf.lattice.num_rings-buffer)]], dims = 2)./(conf.lattice.num_rings-2*buffer)), window_size), label = "13")
for i in 1:13
    lines!(ax, rolling_average(vec(sum(clusters[1,1,:,proto_indices[buffer:(conf.lattice.num_rings-buffer),i]], dims = 2)./(conf.lattice.num_rings-2*buffer)), window_size), label = "$i")

end
axislegend(ax, position=:rc)

display(fig)
#save("clusters.png", fig)

function rolling_average(y, window)
    return [mean(y[max(1, i - window + 1):i]) for i in 1:length(y)]
end
window_size = 1
CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1), xlabel = "timestep", ylabel = "Energy")

lines!(ax, rolling_average(vec(sum(energy_comps[5,9,:,:,seam_left[buffer:(conf.lattice.num_rings-buffer)]], dims = (2,3))./(conf.lattice.num_rings-2*buffer)), window_size), label = "1")
lines!(ax, rolling_average(vec(sum(energy_comps[5,9,:,:,seam_right[buffer:(conf.lattice.num_rings-buffer)]], dims = (2,3))./(conf.lattice.num_rings-2*buffer)), window_size), label = "13")
for i in 2:12
    lines!(ax, rolling_average(vec(sum(energy_comps[5,9,:,:,proto_indices[buffer:(conf.lattice.num_rings-buffer),i]], dims =(2,3))./(conf.lattice.num_rings-2*buffer)), window_size), label = "$i")

end
axislegend(ax, position=:rc)

display(fig)
#save("energy_plot.png", fig)


CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))
hist!(ax,log.(lattice.rates[buffer*conf.lattice.N:(conf.lattice.num_rings-buffer)*conf.lattice.N]))

display(fig)

lattice.rates .> 10

CairoMakie.activate!()

fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1))
hist!(ax,log.(sum(energy_comps[1,1,7201,:,:], dims = 1)[buffer*conf.lattice.N:(conf.lattice.num_rings-buffer)*conf.lattice.N]))

display(fig)


fig = Figure(resolution = (1500, 1000))
ax = Axis(fig[1,1], aspect=AxisAspect(1), xlabel = "Log Transition rate into kinesin-like", ylabel = "Log Transition rate out of kinesin-like", title = "Proportion of tubulin in kinesin-like state at equilibrium")
res = [[mean(vec(sum(clusters[i,j,:,:], dims = 2))[Int64(round(0.9*length(time))):length(time)]) for i in 1:num] for j in 1:num]
#reshape(res, (num,num))
res_new = zeros((num,num))
for i in 1:num
    for j in 1:num
        res_new[i,j] = res[i][j]
    end
end
heatmap!(ax, log10.(a1), log10.(a2), res_new)
Colorbar(fig[:, end+1])
display(fig)
#save("heat_map.png")
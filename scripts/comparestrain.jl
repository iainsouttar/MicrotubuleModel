###one protofilament
protoEnergyArr = zeros(120,13)
protoEnergyArrpref = zeros(120,13)

proto133 = zeros(120,3,13)
for i in 1:120
    for j in 1:13
        #protoEnergyArr[i, j] = sum(E60thirteen3beads[:,:,4,1], dims = 1)[Int((i-1)*13 +j)]
        #protoEnergyArrpref[i, j] = sum(E60thirteen3beads[:,:,4,2], dims = 1)[Int((i-1)*13 +j)]

        proto133[i,:, j] = latticepref120thirteen3x[Int((i-1)*13 +j),:]
    end
end


protoEnergyArr13 = zeros(120,13)
protoEnergyArrpref13 = zeros(120,13)

proto133 = zeros(120,3,13)
proto1332 = zeros(120,3,13)
proto1333 = zeros(120,3,13)
for i in 1:120
    for j in 1:13
        protoEnergyArr13[i, j] = sum(E120thirteenArrbeadslast[:,:], dims = 1)[Int((i-1)*13 +j)]
        protoEnergyArrpref13[i, j] = sum(E120thirteenArrbeadslastpref[:,:], dims = 1)[Int((i-1)*13 +j)]
        proto133[i,:, j] = lattice120thirteenArrxpreflast[Int((i-1)*13 +j),:]
        proto1333[i,:, j] = lattice120thirteenArrx[Int((i-1)*13 +j),:]

        #proto1332[i,:, j] = lattice360thirteen3xlowtime[Int((i-1)*13 +j),:]
    end
end



CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
means = zeros(13)
means1 = zeros(13)
lengths = zeros(13)
lengths1 = zeros(13)
for j in 1:13
    means[j] = mean(protoEnergyArrpref13[40:80,j])
    means1[j] = mean(protoEnergyArr13[40:80,j])
    lengths[j] = norm(proto133[80,:, j] - proto133[40,:,j])
    lengths1[j]=norm(proto1333[80,:, j] - proto1333[40,:,j])

end
lines!(ax, means, label="Preference")
lines!(ax, means1, label="No preference")

ax.xlabel = "Protofilament id"
ax.ylabel = "Mean Energy"

#axislegend(ax, position=:bc)
display(fig)
save("plots/strain13mean.png", fig)

CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])

lines!(ax, lengths, label="Preference")
lines!(ax, lengths1, label="No preference")

ax.xlabel = "Protofilament id"
ax.ylabel = "Length"

#axislegend(ax, position=:bc)
display(fig)
save("plots/strain13lengths.png", fig)




CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
for j in (1,2,3)
    #lines!(ax, protoEnergyArr13[20:100,j], label = "13-3 no pref ID:$j", palette =:tab13)
    lines!(ax, protoEnergyArrpref13[20:100,j], label = "13-3 pref ID:$j", palette =:tab13)

    #lines!(ax, protoEnergyArr14[20:40,j], label = "14-2 no pref ID:$j", palette =:tab13)
    #lines!(ax, protoEnergyArrpref14[20:40,j], label = "14-2 pref ID:$j", palette =:tab13)

end
ax.xlabel = "Monomer number"
ax.ylabel = "Energy"

axislegend(ax, position=:rc)
display(fig)
save("plots/strainpref13.png", fig)


GLMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis3(fig[1,1])
lines!(ax,proto133[2:100,:,1], label = "1", color = :black)
scatter!(ax,proto133[2:2:100,:,1], label = "1", color=:yellow, markersize = 10)
scatter!(ax,proto133[1:2:99,:,1], label = "1", color=:green, markersize = 10)
#lines!(ax,proto1333[2:160,:,1], label = "1", color = :black)
#scatter!(ax,proto1333[2:2:160,:,1], label = "1", color=:yellow, markersize = 10)
#scatter!(ax,proto1333[1:2:159,:,1], label = "1", color=:green, markersize = 10)



for j in 2:13
    lines!(ax,proto133[1:100,:,j], label = "1", color = :black)
    scatter!(ax,proto133[2:2:100,:,j], label = "1", color=:red, markersize = 10)
    scatter!(ax,proto133[1:2:99,:,j], label = "1", color=:blue, markersize = 10)
end

display(fig)



GLMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis3(fig[1,1])
#lines!(ax,proto133[20:220,:,1], label = "1", color = :black)
# scatter!(ax,proto133[20:2:220,:,1], label = "1", color=:yellow, markersize = 10)
# scatter!(ax,proto133[21:2:219,:,1], label = "1", color=:green, markersize = 10)

scatter!(ax,proto133[20:2:100,:,1]-proto1333[20:2:100,:,1], label = "1", color=:yellow, markersize = 10)
scatter!(ax,proto133[19:2:99,:,1]-proto1333[19:2:99,:,1], label = "1", color=:green, markersize = 10)


display(fig)



GLMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis3(fig[1,1])
for j in 1:1
    scatter!(proto133[70:170,:,j])
    scatter!(proto1332[70:170,:,j], color = :red)
    scatter!(proto1333[70:170,:,j], color = :blue)

end

display(fig)



GLMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis3(fig[1,1])

scatter!(ax,lattice60fourteen2x)
scatter!(ax,lattice60fourteen2xlowtime)
display(fig)
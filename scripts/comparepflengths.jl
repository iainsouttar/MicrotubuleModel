#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


Nt = 4_000_000
stp = 1000
time = collect(0:stp:Nt)



conf48lowtime = from_toml(MicrotubuleConfig, "config/eulerMT.toml")
confPref48lowtime = from_toml(MicrotubuleConfig, "config/eulerMTpref.toml")
#conf = set_bond_angles(conf)

lattice48lowtime, bead_info48lowtime= initialise(conf48lowtime)
latticepref48lowtime, bead_infopref48lowtime = initialise(confPref48lowtime)



E48lowtime = zeros((7,length(time),2))
cosangles48lowtime = zeros(length(time), length(latticepref48lowtime),2)
lngth48lowtime = zeros(length(time),2)
E48lowtime[:,1,1], cosangles48lowtime[1, :,1], lngth48lowtime[1,1] = total_energy(lattice48lowtime, bead_info48lowtime)
E48lowtime[:,1,2], cosangles48lowtime[1, :,2], lngth48lowtime[1,2] = total_energy(latticepref48lowtime, bead_infopref48lowtime)


@showprogress for i in 1:Nt
    iterate!(lattice48lowtime, bead_info48lowtime, conf48lowtime, conf48lowtime.iter_pars)
    iterate!(latticepref48lowtime, bead_infopref48lowtime, confPref48lowtime, confPref48lowtime.iter_pars)

    if i % stp == 0
        E48lowtime[:,i÷stp+1,1], cosangles48lowtime[i÷stp+1, :,1], lngth48lowtime[i÷stp+1,1] = total_energy(lattice48lowtime, bead_info48lowtime)
        E48lowtime[:,i÷stp+1,2], cosangles48lowtime[i÷stp+1, :,2], lngth48lowtime[i÷stp+1,2] = total_energy(latticepref48lowtime, bead_infopref48lowtime)

    end
end

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
plot_energy!(ax, time, E48lowtime[:,:,2])
f

save("plots/energyMTprefnoT.png", f)


# CairoMakie.activate!()
# f = Figure(resolution=(1000,600))
# ax = Axis(f[1,1])
# lines!(ax, time, lngth)
# f
# save("plots/MTlengthprefnoT.png", f)

# f = plot_3d(lattice, bead_info)
# f

#Plotting the length of each protofilament making up the mt

distsNopref48lowtime = zeros(conf48lowtime.lattice.N)
distspref48lowtime = zeros(confPref48lowtime.lattice.N)


#3 should be changed if the rise is changed
for i in 1:conf48lowtime.lattice.N
    distsNopref48lowtime[i] = norm(lattice48lowtime.x[i+3*conf48lowtime.lattice.N]-lattice48lowtime.x[(conf48lowtime.lattice.num_rings-4)*conf48lowtime.lattice.N+i])
    distspref48lowtime[i] = norm(latticepref48lowtime.x[i+3*confPref48lowtime.lattice.N]-latticepref48lowtime.x[(confPref48lowtime.lattice.num_rings-4)*confPref48lowtime.lattice.N+i])

end

# for i in 1:conf.lattice.N
#     distsNopref1[i] = norm(lattice.x[i+3*conf.lattice.N]-lattice.x[(conf.lattice.num_rings-4)*conf.lattice.N+i])
#     distspref1[i] = norm(latticepref.x[i+3*confPref.lattice.N]-latticepref.x[(confPref.lattice.num_rings-4)*confPref.lattice.N+i])

# end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf48lowtime.lattice.N, distspref48lowtime, color = :blue)
lines!(ax, 1:conf48lowtime.lattice.N, distsNopref48lowtime, color = :red)

f
save("plots/MTlengthpfs2comp48lowtime.png", f)

distsNopref48lowtime = zeros(conf48lowtime.lattice.N)
distspref48lowtime = zeros(confPref48lowtime.lattice.N)


#3 should be changed if the rise is changed
for i in 1:conf48lowtime.lattice.N
    distsNopref48lowtime[i] = norm(lattice48lowtime.x[i+15*conf48lowtime.lattice.N]-lattice48lowtime.x[(conf48lowtime.lattice.num_rings-16)*conf48lowtime.lattice.N+i])
    distspref48lowtime[i] = norm(latticepref48lowtime.x[i+15*confPref48lowtime.lattice.N]-latticepref48lowtime.x[(confPref48lowtime.lattice.num_rings-16)*confPref48lowtime.lattice.N+i])

end

# for i in 1:conf.lattice.N
#     distsNopref1[i] = norm(lattice.x[i+3*conf.lattice.N]-lattice.x[(conf.lattice.num_rings-4)*conf.lattice.N+i])
#     distspref1[i] = norm(latticepref.x[i+3*confPref.lattice.N]-latticepref.x[(confPref.lattice.num_rings-4)*confPref.lattice.N+i])

# end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf48lowtime.lattice.N, distspref48lowtime, color = :blue)
lines!(ax, 1:conf48lowtime.lattice.N, distsNopref48lowtime, color = :red)

f
save("plots/MTlengthpfs2comp48lowtimecut.png", f)

distsNopref48lowtimewhole = zeros(conf48lowtime.lattice.N)
distspref48lowtimewhole = zeros(confPref48lowtime.lattice.N)
distsdiff48lowtimewhole = zeros(confPref48lowtime.lattice.N)
distsdiff48lowtime = zeros(confPref48lowtime.lattice.N)


#3 should be changed if the rise is changed
for i in 1:conf48lowtime.lattice.N
    distsNopref48lowtimewhole[i] = norm(lattice48lowtime.x[i]-lattice48lowtime.x[(conf48lowtime.lattice.num_rings-1)*conf48lowtime.lattice.N+i])
    distspref48lowtimewhole[i] = norm(latticepref48lowtime.x[i]-latticepref48lowtime.x[(confPref48lowtime.lattice.num_rings-1)*confPref48lowtime.lattice.N+i])
    distsdiff48lowtimewhole[i] = distsNopref48lowtimewhole[i] - distspref48lowtimewhole[i]
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf48lowtime.lattice.N, distspref48lowtimewhole, color = :blue)
lines!(ax, 1:conf48lowtime.lattice.N, distsNopref48lowtimewhole, color = :red)

f
save("plots/MTlengthpfs2comp48lowtimewhole.png", f)




#3 should be changed if the rise is changed
for i in 1:conf48lowtime.lattice.N
    distsdiff48lowtime[i] = distsNopref48lowtime[i]-distspref48lowtime[i]
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf48lowtime.lattice.N, distsdiff48lowtimewhole, color = :blue)

f
save("plots/MTlengthpfs2comp48lowtimewholediff.png", f)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf48lowtime.lattice.N, distsdiff48lowtime, color = :blue)

f
save("plots/MTlengthpfs2comp48lowtimediff.png", f)
#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


Nt = 1000000
stp = 20000
time = collect(0:stp:Nt)



conf120fourteen2 = from_toml(MicrotubuleConfig, "config/eulerMT142.toml")
confPref120fourteen2 = from_toml(MicrotubuleConfig, "config/eulerMTpref142.toml")
#conf = set_bond_angles(conf)

lattice120fourteen2, bead_info120fourteen2= initialise(conf120fourteen2)
latticepref120fourteen2, bead_infopref120fourteen2 = initialise(confPref120fourteen2)



E120fourteen2 = zeros((7,length(time),2))
cosangles120fourteen2 = zeros(length(time), length(latticepref120fourteen2),2)
lngth120fourteen2 = zeros(length(time),2)

E120fourteen2beads = zeros((7,length(lattice120fourteen2), length(time),2))

E120fourteen2[:,1,1], cosangles120fourteen2[1, :,1], lngth120fourteen2[1,1], E120fourteen2beads[:,:,1,1] = total_energy(lattice120fourteen2, bead_info120fourteen2)
E120fourteen2[:,1,2], cosangles120fourteen2[1, :,2], lngth120fourteen2[1,2], E120fourteen2beads[:,:,1,2] = total_energy(latticepref120fourteen2, bead_infopref120fourteen2)


@showprogress for i in 1:Nt
    #iterate!(lattice120fourteen2, bead_info120fourteen2, conf120fourteen2, conf120fourteen2.iter_pars)
    iterate!(latticepref120fourteen2, bead_infopref120fourteen2, confPref120fourteen2, confPref120fourteen2.iter_pars)

    if i % stp == 0
        E120fourteen2[:,i÷stp+1,1], cosangles120fourteen2[i÷stp+1, :,1], lngth120fourteen2[i÷stp+1,1], E120fourteen2beads[:,:,i÷stp+1,1] = total_energy(lattice120fourteen2, bead_info120fourteen2)
        E120fourteen2[:,i÷stp+1,2], cosangles120fourteen2[i÷stp+1, :,2], lngth120fourteen2[i÷stp+1,2], E120fourteen2beads[:,:,i÷stp+1,2] = total_energy(latticepref120fourteen2, bead_infopref120fourteen2)

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
# CairoMakie.activate!()
# f = Figure(resolution=(1000,600))
# ax = Axis(f[1,1])
# plot_energy!(ax, time, E120fourteen2[:,:,2])
# f

# save("plots/energyMTprefnoT.png", f)

writedlm("results/raw/120fourteen2.csv", lattice120fourteen2.x, ',')
writedlm("results/raw/pref120fourteen2.csv", latticepref120fourteen2.x, ',')

writedlm("results/raw/120fourteen2Energy.csv", E120fourteen2beads[:, :,length(time),1], ',')
writedlm("results/raw/120fourteen2EnergyPref.csv", E120fourteen2beads[:,:,length(time),2], ',')


conf240thirteenLArr = from_toml(MicrotubuleConfig, "config/eulerMT.toml")
confPref240thirteenLArr = from_toml(MicrotubuleConfig, "config/eulerMTpref.toml")
#conf = set_bond_angles(conf)

lattice240thirteenLArr, bead_info240thirteenLArr= initialise(conf240thirteenLArr)
latticepref240thirteenLArr, bead_infopref240thirteenLArr = initialise(confPref240thirteenLArr)

E240thirteenLArr = zeros((7,length(time),2))
cosangles240thirteenLArr = zeros(length(time), length(latticepref240thirteenLArr),2)
lngth240thirteenLArr = zeros(length(time),2)

E240thirteenLArrbeads = zeros((7,length(lattice240thirteenLArr), length(time),2))

E240thirteenLArr[:,1,1], cosangles240thirteenLArr[1, :,1], lngth240thirteenLArr[1,1], E240thirteenLArrbeads[:,:,1,1] = total_energy(lattice240thirteenLArr, bead_info240thirteenLArr)
E240thirteenLArr[:,1,2], cosangles240thirteenLArr[1, :,2], lngth240thirteenLArr[1,2], E240thirteenLArrbeads[:,:,1,2] = total_energy(latticepref240thirteenLArr, bead_infopref240thirteenLArr)

lattice240 = []
bead_info240 = []
@showprogress for i in 1:Nt
    #if i%100 ==0
    #    iterate!(lattice240thirteenLArr, bead_info240thirteenLArr, conf240thirteenLArr, conf240thirteenLArr.iter_pars)
    #end
    iterate!(latticepref240thirteenLArr, bead_infopref240thirteenLArr, confPref240thirteenLArr, confPref240thirteenLArr.iter_pars)

    if i % stp == 0
        E240thirteenLArr[:,i÷stp+1,1], cosangles240thirteenLArr[i÷stp+1, :,1], lngth240thirteenLArr[i÷stp+1,1], E240thirteenLArrbeads[:,:,i÷stp+1,1] = total_energy(lattice240thirteenLArr, bead_info240thirteenLArr)
        E240thirteenLArr[:,i÷stp+1,2], cosangles240thirteenLArr[i÷stp+1, :,2], lngth240thirteenLArr[i÷stp+1,2], E240thirteenLArrbeads[:,:,i÷stp+1,2] = total_energy(latticepref240thirteenLArr, bead_infopref240thirteenLArr)
        push!(lattice240, deepcopy(latticepref240thirteenLArr))
        push!(bead_info240, deepcopy(bead_infopref240thirteenLArr))
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
# CairoMakie.activate!()
# f = Figure(resolution=(1000,600))
# ax = Axis(f[1,1])
# plot_energy!(ax, time, E360thirteen3longtime[:,:,2])
# f

# save("plots/energyMTprefnoT.png", f)

writedlm("results/raw/240thirteenLArrEnergy.csv", E240thirteenLArr[:,:,1], ',')
writedlm("results/raw/240thirteenLArrEnergyPref.csv", E240thirteenLArr[:,:,2], ',')


writedlm("results/raw/240thirteenLArrEnergyBeads.csv", E240thirteenLArrbeads[:,:,length(lattice240),1], ',')
writedlm("results/raw/240thirteenLArrEnergyPrefBeads.csv", E240thirteenLArrbeads[:,:,length(lattice240),2], ',')


writedlm("results/raw/240thirteenLArr.csv", lattice240thirteenLArr.x, ',')
writedlm("results/raw/pref240thirteenLArr.csv", latticepref240thirteenLArr.x, ',')


writedlm("results/raw/240thirteenLArr1.csv", lattice240thirteenLArr.x, ',')
writedlm("results/raw/pref240thirteenLArrfirst.csv", lattice240[1].x, ',')
writedlm("results/raw/pref240thirteenLArrlast.csv", lattice240[length(lattice240)].x, ',')
writedlm("results/raw/pref240thirteenLArrSecondlast.csv", lattice240[length(lattice240)-10].x, ',')


lattice240thirteenLArrxpref = readdlm("results/raw/pref240thirteenLArr.csv", ',')
lattice240thirteenLArrx = readdlm("results/raw/240thirteenLArr.csv", ',')

lattice240thirteenLArrxpref1 = readdlm("results/raw/pref240thirteenLArrfirst.csv", ',')
lattice240thirteenLArrxpreflast = readdlm("results/raw/pref240thirteenLArrlast.csv", ',')
lattice240thirteenLArrxprefsecondlast= readdlm("results/raw/pref240thirteenLArrSecondlast.csv", ',')

E240thirteenLArrbeadslast = readdlm("results/raw/240thirteenLArrEnergyBeads.csv", ',')
E240thirteenLArrbeadslastpref =  readdlm("results/raw/240thirteenLArrEnergyPrefBeads.csv", ',')

#lattice240thirteenLArrxpref = readdlm("results/raw/pref240thirteenL3longtimeArr1.csv", ',')



Nt = 3_000000
stp = 100_00
time = collect(0:stp:Nt)



conf24thirteen3kinesin = from_toml(MicrotubuleConfig, "config/eulerMT24.toml")
#conf = set_bond_angles(conf)



lattice24thirteen3kinesin, bead_info24thirteen3kinesin= initialise(conf24thirteen3kinesin)
#latticepref24thirteen3kinesin, bead_infopref24thirteen3kinesin = initialise(confPref24thirteen3kinesin)

attach_kinesin!(lattice24thirteen3kinesin, bead_info24thirteen3kinesin, conf24thirteen3kinesin.lattice.N, conf24thirteen3kinesin.lattice.S, Int(8 + 1))
attach_kinesin!(lattice24thirteen3kinesin, bead_info24thirteen3kinesin, conf24thirteen3kinesin.lattice.N, conf24thirteen3kinesin.lattice.S, Int(16 + 1))


E24thirteen3kinesin = zeros((7,length(time),2))
cosangles24thirteen3kinesin = zeros(length(time), length(lattice24thirteen3kinesin),2)
lngth24thirteen3kinesin = zeros(length(time),2)

E24thirteen3kinesinbeads = zeros((7,length(lattice24thirteen3kinesin), length(time),2))

E24thirteen3kinesin[:,1,1], cosangles24thirteen3kinesin[1, :,1], lngth24thirteen3kinesin[1,1], E24thirteen3kinesinbeads[:,:,1,1] = total_energy(lattice24thirteen3kinesin, bead_info24thirteen3kinesin)
#E24thirteen3kinesin[:,1,2], cosangles24thirteen3kinesin[1, :,2], lngth24thirteen3kinesin[1,2], E24thirteen3kinesinbeads[:,:,1,2] = total_energy(latticepref24thirteen3kinesin, bead_infopref24thirteen3kinesin)


@showprogress for i in 1:Nt
    #iterate!(lattice24thirteen3kinesin, bead_info24thirteen3kinesin, conf24thirteen3kinesin, conf24thirteen3kinesin.iter_pars)
    iterate!(lattice24thirteen3kinesin, bead_info24thirteen3kinesin, conf24thirteen3kinesin, conf24thirteen3kinesin.iter_pars)

    if i % stp == 0
        E24thirteen3kinesin[:,i÷stp+1,1], cosangles24thirteen3kinesin[i÷stp+1, :,1], lngth24thirteen3kinesin[i÷stp+1,1], E24thirteen3kinesinbeads[:,:,i÷stp+1,1] = total_energy(lattice24thirteen3kinesin, bead_info24thirteen3kinesin)
        #E24thirteen3kinesin[:,i÷stp+1,2], cosangles24thirteen3kinesin[i÷stp+1, :,2], lngth24thirteen3kinesin[i÷stp+1,2], E24thirteen3kinesinbeads[:,:,i÷stp+1,2] = total_energy(latticepref24thirteen3kinesin, bead_infopref24thirteen3kinesin)

    end
end


# f = plot_3d(lattice60thirteen3, bead_infopref60thirteen3)
# f

lattice360thirteen3xlowtime = readdlm("results/raw/pref360thirteen3.csv", ',')


lattice60fourteen3x = readdlm("results/raw/60fourteen3.csv", ',')
lattice60twelve3x = readdlm("results/raw/60twelve3.csv", ',')
latticepref60fourteen3x = readdlm("results/raw/pref60fourteen3.csv", ',')
latticepref60twelve3x = readdlm("results/raw/pref60twelve3.csv", ',')

lattice60thirteen3x = readdlm("results/raw/60thirteen3.csv", ',')
lattice60fourteen2x = readdlm("results/raw/60fourteen2.csv", ',')
latticepref60thirteen3x = readdlm("results/raw/pref60thirteen3.csv", ',')
latticepref60fourteen2x = readdlm("results/raw/pref60fourteen2.csv", ',')


latticepref240thirteenL3x = readdlm("results/raw/pref240thirteenL3.csv", ',')


midPointRings143 = zeros(60, 3)
midPointRings143Pref = zeros(60, 3)
midPointRingsPref123 =zeros(60,3)
midPointRings123 = zeros(60,3)
for i in 1:60
    midPointRings123[i,:] = sum(lattice60twelve3x[Int((i-1)*12 +1):Int((i)*12),:], dims = 1)/12 
    midPointRingsPref123[i,:] = sum(latticepref60twelve3x[Int((i-1)*12 +1):Int((i)*12),:], dims = 1)/12 


    midPointRings142[i,:] = sum(lattice60fourteen3x[Int((i-1)*14 +1):Int((i)*14),:], dims = 1)/14 
    midPointRings142Pref[i,:] = sum(latticepref60fourteen3x[Int((i-1)*14 +1):Int((i)*14),:], dims = 1)/14 

end


using GLMakie
# Create figure and axes
fig = Figure(resolution = (1000, 1500))
ax = Axis3(fig[1,1])

# Create scatter plot
scatter!(ax,midPointRings123[30:40, :], color = :blue)
scatter!(ax, midPointRings143[30:40, :], color = :red)


# Set title
#fig[1, 1] = title("3D Scatter Plot")

# Show the plot
display(fig)


diffs123 = zeros(Int(60-1),3)
diffs133 = zeros(Int(60-1),3)
for i in 1:Int(60-1)
    diffs123[i,:] = midPointRings123[i+1,:]-midPointRings123[i,:]
    diffs133[i,:] = midPointRings133[i+1,:]-midPointRings133[i,:]
end


###one protofilament
proto142 = zeros(60,3,13)
proto2142 = zeros(60,3)
skew = zeros()
for i in 1:60
    for j in 1:13
        proto142[i,:, j] = lattice60fourteen2x[Int((i-1)*13 +j),:]
    end
end




fig142 = Figure(resolution = (1000, 750))
ax = Axis3(fig142[1,1])
for j in 1:13
    scatter!(ax,proto142[:,:,j], label = "$j")
end
# axislegend(ax, position=:lc)
# ax.xlabel = "Longitudinal id"
# ax.ylabel = "Angle between bond below and bond above (radians)"
# ax.title = "N = 14, S = 2, without preference"
display(fig142)
#save("plots/protofilament60fourteen2.png", fig)


display(fig143)
display(fig133)
display(fig142)

display(fig133pref)



diffs142 = zeros(Int(60-1),3,14)
angles = zeros(Int(60-2),14)
for i in 1:Int(60-1)
    for j in 1:14
        diffs142[i,:,j] = proto142[i+1,:,j] - proto142[i,:,j]
        if i>1
            angles[i-1,j] = (acos(dot(diffs142[i,:,j], diffs142[i-1,:,j])/(norm(diffs143[i,:,j])*norm(diffs143[i-1,:,j]))))
        end
    end
end



fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
for j in 1:14
    lines!(ax, angles[20:40,j], label = "$j", palette =:tab10)
end
axislegend(ax, position=:lc)
ax.xlabel = "Longitudinal id"
ax.ylabel = "Angle between bond below and bond above (radians)"
ax.title = "N = 14, S = 2, without preference"
display(fig)
save("plots/Skew_angle_60fourteen2.png", fig)

using GLMakie
# Create figure and axes
fig = Figure(resolution = (1000, 1500))
ax = Axis3(fig[1,1])

# Create scatter plot
scatter!(ax,diffs123, color = :blue)


# Set title
#fig[1, 1] = title("3D Scatter Plot")

# Show the plot
display(fig)


# for i in 2:60
#     #calc skew
#     skew[i] = proto[i-1]
# end
fig = Figure(resolution = (1000, 1500))
ax = Axis3(fig[1,1])

# Create scatter plot
scatter!(ax,proto[30:40], color = :blue)


# Set title
#fig[1, 1] = title("3D Scatter Plot")

# Show the plot
display(fig)



#scatter(diffs142[20:30,:])

# f = plot_3d(lattice60fourteen3, bead_info60fourteen3)
# f

using GLMakie
# Create figure and axes
fig = Figure(resolution = (1000, 1500))
ax = Axis3(fig[1,1])

# Create scatter plot
scatter!(ax,lattice60fourteen2x, color = :blue)
#scatter!(ax, midPointRings142, color = :red)


# Set title
#fig[1, 1] = title("3D Scatter Plot")

# Show the plot
display(fig)


distsNopref60twelve3 = zeros(conf60twelve3.lattice.N)
distspref60twelve3 = zeros(confPref60twelve3.lattice.N)


#3 should be changed if the rise is changed
for i in 1:conf60twelve3.lattice.N
    distsNopref60twelve3[i] = norm(lattice60twelve3.x[i+25*conf60twelve3.lattice.N]-lattice60twelve3.x[(conf60twelve3.lattice.num_rings-25)*conf60twelve3.lattice.N+i])
    distspref60twelve3[i] = norm(latticepref60twelve3.x[i+25*confPref60twelve3.lattice.N]-latticepref60twelve3.x[(confPref60twelve3.lattice.num_rings-25)*confPref60twelve3.lattice.N+i])

end

# for i in 1:conf.lattice.N
#     distsNopref1[i] = norm(lattice.x[i+3*conf.lattice.N]-lattice.x[(conf.lattice.num_rings-4)*conf.lattice.N+i])
#     distspref1[i] = norm(latticepref.x[i+3*confPref.lattice.N]-latticepref.x[(confPref.lattice.num_rings-4)*confPref.lattice.N+i])

# end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf60twelve3.lattice.N, distspref60twelve3, color = :blue, label = "Difference in preference")
lines!(ax, 1:conf60twelve3.lattice.N, distsNopref60twelve3, color = :red, label = "No difference")
ax.ylabel = "Length"
ax.xlabel = "Protofilament id"
axislegend(ax, position=:cc)


f
save("plots/MTlengthpfs2comp60twelve3cut.png", f)


Nt = 6_00000
stp = 100_00
timex = collect(0:stp:Nt)



conf24thirteen3kinesinReal = from_toml(MicrotubuleConfig, "config/eulerMTpref.toml")
#conf = set_bond_angles(conf)



lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal= initialise(conf24thirteen3kinesinReal)
#latticepref24thirteen3kinesinReal, bead_infopref24thirteen3kinesinReal = initialise(confPref24thirteen3kinesinReal)

attach_kinesin!(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal, conf24thirteen3kinesinReal.lattice.N, conf24thirteen3kinesinReal.lattice.S, Int(130 + 5))
#attach_kinesin!(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal, conf24thirteen3kinesinReal.lattice.N, conf24thirteen3kinesinReal.lattice.S, Int( + 1))


E24thirteen3kinesinReal = zeros((7,length(timex),2))
cosangles24thirteen3kinesinReal = zeros(length(timex), length(lattice24thirteen3kinesinReal),2)
lngth24thirteen3kinesinReal = zeros(length(timex),2)

E24thirteen3kinesinRealbeads = zeros((7,length(lattice24thirteen3kinesinReal), length(timex),2))

E24thirteen3kinesinReal[:,1,1], cosangles24thirteen3kinesinReal[1, :,1], lngth24thirteen3kinesinReal[1,1], E24thirteen3kinesinRealbeads[:,:,1,1] = total_energy(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal)
#E24thirteen3kinesinReal[:,1,2], cosangles24thirteen3kinesinReal[1, :,2], lngth24thirteen3kinesinReal[1,2], E24thirteen3kinesinRealbeads[:,:,1,2] = total_energy(latticepref24thirteen3kinesinReal, bead_infopref24thirteen3kinesinReal)


@showprogress for i in 1:Nt
    #iterate!(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal, conf24thirteen3kinesinReal, conf24thirteen3kinesinReal.iter_pars)
    iterate!(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal, conf24thirteen3kinesinReal, conf24thirteen3kinesinReal.iter_pars)

    if i % stp == 0
        E24thirteen3kinesinReal[:,i÷stp+1,1], cosangles24thirteen3kinesinReal[i÷stp+1, :,1], lngth24thirteen3kinesinReal[i÷stp+1,1], E24thirteen3kinesinRealbeads[:,:,i÷stp+1,1] = total_energy(lattice24thirteen3kinesinReal, bead_info24thirteen3kinesinReal)
        #E24thirteen3kinesinReal[:,i÷stp+1,2], cosangles24thirteen3kinesinReal[i÷stp+1, :,2], lngth24thirteen3kinesinReal[i÷stp+1,2], E24thirteen3kinesinRealbeads[:,:,i÷stp+1,2] = total_energy(latticepref24thirteen3kinesinReal, bead_infopref24thirteen3kinesinReal)

    end
end


writedlm("results/raw/24thirteen3kinesinRealEnergy.csv", E24thirteen3kinesinReal[:,:,1], ',')


writedlm("results/raw/24thirteen3kinesinRealEnergyBeads.csv", E24thirteen3kinesinRealbeads[:,:,length(timex),1], ',')


writedlm("results/raw/24thirteen3kinesinReal.csv", lattice24thirteen3kinesinReal.x, ',')




lattice24thirteen3kinesinRealx = readdlm("results/raw/24thirteen3kinesinReal.csv", ',')


E24thirteen3kinesinRealbeadslast = readdlm("results/raw/24thirteen3kinesinRealEnergyBeads.csv", ',')



conf24thirteen3 = from_toml(MicrotubuleConfig, "config/eulerMTpref.toml")
#conf = set_bond_angles(conf)



lattice24thirteen3, bead_info24thirteen3= initialise(conf24thirteen3)
#latticepref24thirteen3, bead_infopref24thirteen3 = initialise(confPref24thirteen3)

#attach_kinesin!(lattice24thirteen3, bead_info24thirteen3, conf24thirteen3.lattice.N, conf24thirteen3.lattice.S, Int(130 + 5))
#attach_kinesin!(lattice24thirteen3, bead_info24thirteen3, conf24thirteen3.lattice.N, conf24thirteen3.lattice.S, Int( + 1))


E24thirteen3 = zeros((7,length(timex),2))
cosangles24thirteen3 = zeros(length(timex), length(lattice24thirteen3),2)
lngth24thirteen3 = zeros(length(timex),2)

E24thirteen3beads = zeros((7,length(lattice24thirteen3), length(timex),2))

E24thirteen3[:,1,1], cosangles24thirteen3[1, :,1], lngth24thirteen3[1,1], E24thirteen3beads[:,:,1,1] = total_energy(lattice24thirteen3, bead_info24thirteen3)
#E24thirteen3[:,1,2], cosangles24thirteen3[1, :,2], lngth24thirteen3[1,2], E24thirteen3beads[:,:,1,2] = total_energy(latticepref24thirteen3, bead_infopref24thirteen3)


@showprogress for i in 1:Nt
    #iterate!(lattice24thirteen3, bead_info24thirteen3, conf24thirteen3, conf24thirteen3.iter_pars)
    iterate!(lattice24thirteen3, bead_info24thirteen3, conf24thirteen3, conf24thirteen3.iter_pars)

    if i % stp == 0
        E24thirteen3[:,i÷stp+1,1], cosangles24thirteen3[i÷stp+1, :,1], lngth24thirteen3[i÷stp+1,1], E24thirteen3beads[:,:,i÷stp+1,1] = total_energy(lattice24thirteen3, bead_info24thirteen3)
        #E24thirteen3[:,i÷stp+1,2], cosangles24thirteen3[i÷stp+1, :,2], lngth24thirteen3[i÷stp+1,2], E24thirteen3beads[:,:,i÷stp+1,2] = total_energy(latticepref24thirteen3, bead_infopref24thirteen3)

    end
end


writedlm("results/raw/24thirteen3Energy.csv", E24thirteen3[:,:,1], ',')


writedlm("results/raw/24thirteen3EnergyBeads.csv", E24thirteen3beads[:,:,length(timex),1], ',')


writedlm("results/raw/24thirteen3.csv", lattice24thirteen3.x, ',')




lattice24thirteen3x = readdlm("results/raw/24thirteen3.csv", ',')


E24thirteen3beadslast = readdlm("results/raw/24thirteen3EnergyBeads.csv", ',')





conf24thirteen3kinesinSeam = from_toml(MicrotubuleConfig, "config/eulerMTpref.toml")
#conf = set_bond_angles(conf)



lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam= initialise(conf24thirteen3kinesinSeam)
#latticepref24thirteen3kinesinSeam, bead_infopref24thirteen3kinesinSeam = initialise(confPref24thirteen3kinesinSeam)

attach_kinesin!(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam, conf24thirteen3kinesinSeam.lattice.N, conf24thirteen3kinesinSeam.lattice.S, Int(130 + 5))
#attach_kinesin!(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam, conf24thirteen3kinesinSeam.lattice.N, conf24thirteen3kinesinSeam.lattice.S, Int( + 1))


E24thirteen3kinesinSeam = zeros((7,length(timex),2))
cosangles24thirteen3kinesinSeam = zeros(length(timex), length(lattice24thirteen3kinesinSeam),2)
lngth24thirteen3kinesinSeam = zeros(length(timex),2)

E24thirteen3kinesinSeambeads = zeros((7,length(lattice24thirteen3kinesinSeam), length(timex),2))

E24thirteen3kinesinSeam[:,1,1], cosangles24thirteen3kinesinSeam[1, :,1], lngth24thirteen3kinesinSeam[1,1], E24thirteen3kinesinSeambeads[:,:,1,1] = total_energy(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam)
#E24thirteen3kinesinSeam[:,1,2], cosangles24thirteen3kinesinSeam[1, :,2], lngth24thirteen3kinesinSeam[1,2], E24thirteen3kinesinSeambeads[:,:,1,2] = total_energy(latticepref24thirteen3kinesinSeam, bead_infopref24thirteen3kinesinSeam)


@showprogress for i in 1:Nt
    #iterate!(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam, conf24thirteen3kinesinSeam, conf24thirteen3kinesinSeam.iter_pars)
    iterate!(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam, conf24thirteen3kinesinSeam, conf24thirteen3kinesinSeam.iter_pars)

    if i % stp == 0
        E24thirteen3kinesinSeam[:,i÷stp+1,1], cosangles24thirteen3kinesinSeam[i÷stp+1, :,1], lngth24thirteen3kinesinSeam[i÷stp+1,1], E24thirteen3kinesinSeambeads[:,:,i÷stp+1,1] = total_energy(lattice24thirteen3kinesinSeam, bead_info24thirteen3kinesinSeam)
        #E24thirteen3kinesinSeam[:,i÷stp+1,2], cosangles24thirteen3kinesinSeam[i÷stp+1, :,2], lngth24thirteen3kinesinSeam[i÷stp+1,2], E24thirteen3kinesinSeambeads[:,:,i÷stp+1,2] = total_energy(latticepref24thirteen3kinesinSeam, bead_infopref24thirteen3kinesinSeam)

    end
end


writedlm("results/raw/24thirteen3kinesinSeamEnergy.csv", E24thirteen3kinesinSeam[:,:,1], ',')


writedlm("results/raw/24thirteen3kinesinSeamEnergyBeads.csv", E24thirteen3kinesinSeambeads[:,:,length(timex),1], ',')


writedlm("results/raw/24thirteen3kinesinSeam.csv", lattice24thirteen3kinesinSeam.x, ',')




lattice24thirteen3kinesinSeamx = readdlm("results/raw/24thirteen3kinesinSeam.csv", ',')


E24thirteen3kinesinSeambeadslast = readdlm("results/raw/24thirteen3kinesinSeamEnergyBeads.csv", ',')





proto133 = zeros(120,3,13)
proto1332 = zeros(120,3,13)
proto1333 = zeros(120,3,13)
protoEnergyArr13 = zeros(24,13)

for i in 1:24
    for j in 1:13
        protoEnergyArr13[i, j] = sum(E24thirteen3beads[:,:, length(timex),1], dims = 1)[Int((i-1)*13 +j)]
        #protoEnergyArrpref13[i, j] = sum(E24thirteen3kinesinbeadspref[:,:, length(time), 1], dims = 1)[Int((i-1)*13 +j)]
        #proto133[i,:, j] = lattice24thirteen3kinesinxpreflast[Int((i-1)*13 +j),:]
        proto1333[i,:, j] = lattice24thirteen3kinesinRealx[Int((i-1)*13 +j),:]

        #proto1332[i,:, j] = lattice360thirteen3xlowtime[Int((i-1)*13 +j),:]
    end
end

CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1], limits = ((8,16), nothing))
for j in (1,2,13,12, 3,11)
    #lines!(ax, protoEnergyArr13[20:100,j], label = "13-3 no pref ID:$j", palette =:tab13)
    lines!(ax, protoEnergyArr13[:,j], label = "13-3 pref ID:$j", palette =:tab13)

    #lines!(ax, protoEnergyArr14[20:40,j], label = "14-2 no pref ID:$j", palette =:tab13)
    #lines!(ax, protoEnergyArrpref14[20:40,j], label = "14-2 pref ID:$j", palette =:tab13)

end
ax.xlabel = "Monomer number"
ax.ylabel = "Energy"


axislegend(ax, position=:rc)

display(fig)

using GLMakie
GLMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis3(fig[1,1])
for j in (5)
    #lines!(ax, protoEnergyArr13[20:100,j], label = "13-3 no pref ID:$j", palette =:tab13)
    #lines!(ax, protoEnergyArr13[:,j], label = "13-3 pref ID:$j", palette =:tab13)
    lines!(proto1333[8:17, :,j], label = "ID:$j")
    #lines!(ax, protoEnergyArr14[20:40,j], label = "14-2 no pref ID:$j", palette =:tab13)
    #lines!(ax, protoEnergyArrpref14[20:40,j], label = "14-2 pref ID:$j", palette =:tab13)

end


#axislegend(ax, position=:rc)
display(fig)


using GLMakie
# Create figure and axes
# fig = Figure(resolution = (1000, 1500))
# ax = Axis3(fig[1,1])

# # Create scatter plot
# scatter!(ax,lattice24thirteen3kinesin.x)
# #scatter!(ax, midPointRings142, color = :red)


# # Set title
# #fig[1, 1] = title("3D Scatter Plot")

# # Show the plot
# display(fig)

f = plot_3d(lattice24thirteen3kinesin, bead_info24thirteen3kinesin)
f



CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
means = zeros(13)
means1 = zeros(13)
lengths = zeros(13)
lengths1 = zeros(13)
for j in (4,5,6)
    #means[j] = mean(protoEnergyArrpref13[8:16,j])
    means1[j] = mean(protoEnergyArr13[8:12,j])
    #lengths[j] = norm(proto133[80,:, j] - proto133[40,:,j])
    lengths1[j]=norm(proto1333[16,:, j] - proto1333[8,:,j])
end
#lines!(ax, means, label="Preference")
lines!(ax, lengths1[4:6], label="13")

ax.xlabel = "Protofilament id"
ax.ylabel = "Mean Energy"

#axislegend(ax, position=:bc)
display(fig)
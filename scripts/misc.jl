# #Plotting the length of each protofilament making up the mt
# distsNopref1 = zeros(conf.lattice.N)
# for i in 1:conf.lattice.N
#     distsNopref1[i] = norm(lattice.x[i]-lattice.x[(conf.lattice.num_rings-1)*conf.lattice.N+i])
# end

# CairoMakie.activate!()
# f = Figure(resolution=(1000,600))
# ax = Axis(f[1,1])
# lines!(ax, 1:conf.lattice.N, distspref1)
# lines!(ax, 1:conf.lattice.N, distsNopref1)

# f
# save("plots/MTlengthpfs2comp.png", f)


#Plotting the length of each protofilament making up the mt

distsNopref24whole = zeros(conf24.lattice.N)
distspref24whole = zeros(confPref24.lattice.N)
distsdiff24whole = zeros(confPref24.lattice.N)
distsdiff24 = zeros(confPref24.lattice.N)


#3 should be changed if the rise is changed
for i in 1:conf24.lattice.N
    distsNopref24whole[i] = norm(lattice24.x[i]-lattice24.x[(conf24.lattice.num_rings-1)*conf24.lattice.N+i])
    distspref24whole[i] = norm(latticepref24.x[i]-latticepref24.x[(confPref24.lattice.num_rings-1)*confPref24.lattice.N+i])
    distsdiff24whole[i] = distsNopref24whole[i] - distspref24whole[i]
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf24.lattice.N, distspref24whole, color = :blue)
lines!(ax, 1:conf24.lattice.N, distsNopref24whole, color = :red)

f
save("plots/MTlengthpfs2comp24whole.png", f)




#3 should be changed if the rise is changed
for i in 1:conf24.lattice.N
    distsdiff24[i] = distsNopref24[i]-distspref24[i]
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf24.lattice.N, distsdiff24whole, color = :blue)

f
save("plots/MTlengthpfs2comp24wholediff.png", f)

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 1:conf24.lattice.N, distsdiff24, color = :blue)

f
save("plots/MTlengthpfs2comp24diff.png", f)


#plotting protofilaments
protospref2414whole = zeros(conf2414.lattice.num_rings-3,3)
protos2414whole = zeros(conf2414.lattice.num_rings-3,3)

for i in 3:(conf2414.lattice.num_rings-1)
    protos2414whole[i-2,:] = lattice2414.x[(i)*(conf2414.lattice.N)+1]
    protospref2414whole[i-2,:] = latticepref2414.x[(i)*(conf2414.lattice.N)+1]


end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, protospref2414whole[:, 1], protospref2414whole[:,2], color = :blue)
lines!(ax, protos2414whole[:,1], protos2414whole[:,2], color = :red)

f
save("plots/MTlengthpfs2comp2414wholeprotosxy.png", f)
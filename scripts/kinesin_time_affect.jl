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
#noise_list = [0.0, 0.04,0.004]
noise_list = [0.0004, 0.004]
factors = [0.1,1.0,10.0]
rate_list = [0.0,0.0000000000001,0.000000000001, 0.00000000001, 0.0000000001]

noise_list = [0.004]
factors = [0.1, 0.01]
rate_list = [0.0, 0.0000001, 0.000000001,0.00000000001, 0.0000000000001]

noise_list = [0.004]
factors = [1e-8]
rate_list = [1e-8] #, 1e-14]
#0 single ring, 1 big ring, 2 protofilaments, 3 patch, 4 single
kin_shape_list = [4]
kin_start = 5000
kin_taken = 40000


lengthsreps = zeros(reps, length(time), length(noise_list), length(rate_list), length(factors), length(kin_shape_list))
positions = []
lattice, bead_info = initialise(conf)
kinesins = falses(length(time) ,length(lattice.x), reps, length(noise_list), length(rate_list), length(factors), length(kin_shape_list))
clusters = falses(length(time) ,length(lattice.x), reps, length(noise_list), length(rate_list), length(factors), length(kin_shape_list))
lattice_rates = zeros(conf.lattice.num_rings*13, length(time), reps, length(noise_list), length(rate_list), length(factors), length(kin_shape_list))
energies  = zeros(conf.lattice.num_rings*13, length(time), reps, length(noise_list), length(rate_list), length(factors), length(kin_shape_list))

#print(lattice_rates)
for (p,kin_shape) in enumerate(kin_shape_list)
    for (n, factor) in enumerate(factors)
        for (m, trans_rates) in enumerate(rate_list)
            for (l, noises) in enumerate(noise_list)
                for rep in 1:Int(reps)

                    #lattice, bead_info = initialise(conf)
                    #print(lattice.rates)

                    # for i in 1:13
                    #     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(10*13+i), true)
                    # end
                    #positions = []
                    # for j in 1:13
                    #     attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(10*13+j), true)
                    #     print(10*13+j, lattice.kinesin[10*13+j])
                    # end

                    lengths_protos = zeros(length(time), conf.lattice.N)
                    for j in 1:conf.lattice.N
                        lengths_protos[1, j] = norm(lattice.x[5*conf.lattice.N + j]-lattice.x[(conf.lattice.num_rings-5)*conf.lattice.N + j])
                    end
                    @showprogress for i in 1:Nt
                        #iterate!(lattice120fourteen2, bead_info120fourteen2, conf120fourteen2, conf120fourteen2.iter_pars)
                        iterate!(lattice, bead_info, conf, conf.iter_pars,0, trans_rates, trans_rates*factor, noises)
                        #0.00000000001
                        if i % stp == 0
                            kinesins[Int(i/stp) + 1,:,rep, l,m,n,p] = lattice.kinesin[:,1]
                            clusters[Int(i/stp) + 1,:,rep,l,m,n,p] = lattice.kinesin[:,2]
                            for j in 1:conf.lattice.N
                                lengths_protos[Int(i/stp) + 1, j] = norm(lattice.x[10*conf.lattice.N+j]-lattice.x[(conf.lattice.num_rings-10)*conf.lattice.N + j])
                            end
                            lattice_rates[:,Int(i/stp) + 1, rep, l,m,n,p] = lattice.rates[:]
                            #energies[:,Int(i/stp) + 1, rep, l,m,n,p] = 
                            push!(positions, deepcopy(lattice.x))
                        end
                        if i == kin_start
                            if kin_shape == 0
                                for k in 0:0
                                    for j in 1:13
                                        attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false)
                                    end
                                end
                            elseif kin_shape == 1
                                for k in -3:3
                                    for j in 1:13
                                        attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false)
                                    end
                                end
                            elseif kin_shape == 2
                                for k in -(conf.lattice.num_rings/2 -3):(conf.lattice.num_rings/2 -3)
                                    for j in 1:1
                                        attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false)
                                    end
                                end

                            elseif kin_shape == 3
                                for k in -3:3
                                    for j in 1:3
                                        attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false)
                                    end
                                end 
                            else 
                                for k in 0:0
                                    for j in 4:4
                                        attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false)
                                    end
                                end 
                            end
                        end
                        if i == kin_taken
                            if kin_shape == 0
                                for k in 0:0
                                    for j in 1:13
                                        remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false, trans_rates==0.0)
                                    end
                                end
                            elseif kin_shape == 1
                                for k in -3:3
                                    for j in 1:13
                                        remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false, trans_rates==0.0)
                                    end
                                end
                            elseif kin_shape == 2
                                for k in -(conf.lattice.num_rings/2 -3):(conf.lattice.num_rings/2 -3)
                                    for j in 1:1
                                        remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false, trans_rates==0.0)
                                    end
                                end

                            elseif kin_shape == 3
                                for k in -3:3
                                    for j in 1:3
                                        remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false, trans_rates==0.0)
                                    end
                                end 
                            else 
                                for k in 0:0
                                    for j in 4:4
                                        remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((conf.lattice.num_rings/2+k)*13+j), false, trans_rates==0.0)
                                    end
                                end 
                            end
                        end


                    end
                    lengthsreps[rep,:,l,m,n,p] = mean(eachcol(lengths_protos))
                end
            end
        end
    end
end
# using JLD
# save("kinesins.jld","kinesins", kinesins)
# save("clusters.jld","clusters", clusters)
# save("lengthsreps.jld","lengthsreps", lengthsreps)



CairoMakie.activate!()


fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
kin_shapes = ["Single ring", "Large ring", "Protofilaments", "Patch"]
#lines!(ax, mean(eachrow(lengthsreps[:,2:length(time),2])), label = "low noise")
for i in 1:4
    for j in 1:1
        # fig = Figure(resolution = (1000, 750))
        # ax = Axis(fig[1,1])
        #rate_fixed = rate_list[i]
        shape_fixed = kin_shapes[i]
        #factor_fixed = factors[j]
        data = mean(eachcol(vec(sum(clusters[16000:30000,:,1,1, 1, j, i], dims = 2))))
        if i==2
            #data[6000:8000]=repeat([6], 2001)
            #data[8000:9500]=repeat([4], 1501)
            #data[9500:14000] = repeat([0], 4501)
        end
        if i==3
            #data = repeat(data,2)
            #data[6000:8000]=repeat([6], 2001)
            #data[8000:9500]=repeat([4], 1501)
            #data[9500:14000] = repeat([0], 4501)
            #data = repeat(data,inner=[2])

        end
        lines!(ax, collect(0:stp*0.2/(1e2):(stp*0.2*(14001-2))/(1e2)), data[1:14000]./data[1], label = "$shape_fixed")
#        display(fig)
    end
end        
axislegend(ax, position=:rt)
ax.xlabel = "Time after washing out kinesin (microseconds)" 
ax.ylabel = "Proportion of tubulin in kinesin-like conformation"
save("number_decay_shapes.png",fig)
display(fig)
#clusters_low = clusters

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
rate_fixed1 = rate_list[2]
rate_fixed2 = rate_list[3]
lbls = ["Two state", "Three state, $rate_fixed1", "Three state, $rate_fixed2"]
#lines!(ax, mean(eachrow(lengthsreps[:,2:length(time),2])), label = "low noise")
for i in 1:3
    for j in 1:1
        
        rate_fixed = rate_list[i]
        #factor_fixed = factors[j]
        data =  mean(eachrow(lengthsreps[1:reps,2:length(time),1, i, j,2]))
        if i == 3
            data[15999:(length(time)-1)] = repeat(data[15999:(length(time)-1)],inner=[3])[1:Int((length(time) - 15999))]
        end
        lines!(ax, collect(0:stp*0.2/(1e2):(stp*0.2*(length(time)-2))/(1e2)),data[1:(length(time)-1)], label = lbls[i])
    end
end

strng_add = ""
#lines!([collect(0:16/1000:(16*(length(time)-2))/1000)[500],collect(0:16/1000:(16*(length(time)-2))/1000)[500]], [113,115], label = "Kinesin added")
#lines!([collect(0:16/1000:(16*(length(time)-2))/1000)[2500],collect(0:16/1000:(16*(length(time)-2))/1000)[2500]], [113,115], label = "Kinesin added")


    #lines!(ax, collect(0:16/1000:(16*899)/1000), mean(eachrow(lengthsreps[1:reps,2:length(time),1])), label = "high noise")

axislegend(ax, position=:rt)
ax.xlabel = "Time after kinesin binding (microseconds)"
ax.ylabel = "Length of protofilament section (nm)"
ax.title = "Lengths at different rates"
limits!(5,700,113.25,114.3)
display(fig)
save(string("plots/length_relaxation_",strng_add,".png"), fig)

CairoMakie.activate!()

for j in 1:1
    fig = Figure(resolution = (1000, 750))
    ax = Axis(fig[1,1])
    # plot!(ax, clusters[:,12*13:15*13,1,1,2,1,2])
    plot!(ax, clusters[:,:,1,1,j,1,1])
    display(fig)
end
save(string("plots/switching_plot_",strng_add,".png"), fig)

CairoMakie.activate!()
for j in 1:1
    fig = Figure(resolution = (1000, 750))
    ax = Axis(fig[1,1])
    plot!(ax, (log10.((1/rate_list[1])*transpose((lattice_rates[300:350,1:stp:(length(time)),1,1,j,1,1])))))
    #ax.title()
    display(fig)
end
save(string("plots/switching_rates_",strng_add,".png"), fig)

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
for j in 0:4
    for k in 25:30
        lines!(ax,collect(0:stp*0.2/(1e3):(stp*0.2*(length(time)-2))/(1e3))[100:3000], log10.((1/rate_list[1])*lattice_rates[k*13+4+j,2:length(time),1,1,1,1,1])[100:3000])
    end
end
lines!(ax,collect(0:stp*0.2/(1e3):(stp*0.2*(length(time)-2))/(1e3))[100:3000], log10.((1/rate_list[1])*lattice_rates[25*13+4+0,2:length(time),1,1,1,1,1])[100:3000], label = "kinesin bound tubulin")

ax.xlabel ="Time (microseconds)"
#lines!(ax, log10.((1/rate_list[1])*lattice_rates[8*13+4,2:length(time),1,1,1,1,1]))
axislegend(ax, position=:rc)
ax.ylabel = "Log10(Energy)"
ax.title = "Spread of strain"
save("together_inesin_strain.png", fig)
display(fig)



function plot_3d(lattice, bead_info)
    GLMakie.activate!()
    GLMakie.closeall()
    scene = plot(lattice, bead_info)
    return scene
end

f = plot_3d(lattice, bead_info)
f

lines([positions[450][conf.lattice.N*i][1] for i in 5:19])

#lengthsreps = lengthsreps_noswitch
neighbours_indices = zeros(length(time), conf.lattice.num_rings*13)
for t in 1:length(time)
    for i in 1:conf.lattice.num_rings*13
        for bond in bead_info[i].bonds
            if kinesins[t,bond,1,1,1,1,1]
                neighbours_indices[t,i] = 2
                if !kinesins[t,i,1,1,1,1,1]
                    neighbours_indices[t,i] = 1
                end
            end
        end
    end
end




hist(log10.((1/rate_list[1])lattice_rates[neighbours_indices[40000,:].==0,40000,1,1,1,1,1]))
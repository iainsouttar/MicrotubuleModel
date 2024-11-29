if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


Nt = 120000
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
factors = [0.0]
rate_list = [0.0] #, 1e-14]


lengthsreps = zeros(reps, length(time), length(noise_list), length(rate_list), length(factors))
positions = []
lattice, bead_info = initialise(conf)
kinesins = falses(length(time) ,length(lattice.x), reps, length(noise_list), length(rate_list), length(factors))
clusters = falses(length(time) ,length(lattice.x), reps, length(noise_list), length(rate_list), length(factors))
lattice_rates = zeros(48*13, length(time), reps, length(noise_list), length(rate_list), length(factors))
#print(lattice_rates)

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
                    lengths_protos[1, j] = norm(lattice.x[10*conf.lattice.N + j]-lattice.x[(conf.lattice.num_rings-10)*conf.lattice.N + j])
                end
                @showprogress for i in 1:Nt
                    #iterate!(lattice120fourteen2, bead_info120fourteen2, conf120fourteen2, conf120fourteen2.iter_pars)
                    iterate!(lattice, bead_info, conf, conf.iter_pars,0, trans_rates, trans_rates*factor, noises)
                    #0.00000000001
                    if i % stp == 0
                        kinesins[Int(i/stp) + 1,:,rep, l,m,n] = lattice.kinesin[:,1]
                        clusters[Int(i/stp) + 1,:,rep,l,m,n] = lattice.kinesin[:,2]
                        for j in 1:conf.lattice.N
                            lengths_protos[Int(i/stp) + 1, j] = norm(lattice.x[10*conf.lattice.N+j]-lattice.x[(conf.lattice.num_rings-10)*conf.lattice.N + j])
                        end
                        lattice_rates[:,Int(i/stp) + 1, rep, l,m,n] = lattice.rates[:]
                        push!(positions, deepcopy(lattice.x))
                    end
                    for k in -3:3
                        for j in 1:13
                            if i == 2000 #+ j*200
                                    #print(10*13+j, lattice.kinesin[10*13+j])
                                attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((24+k)*13+j), false)
                                    #attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(22*13+j), false)
                            end
                        end
                    end
                    for k in -3:3
                        for j in 1:13
                            if i == 40000 #-j*100
                                #print(10*13+j, lattice.kinesin[10*13+j])
                                remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int((24+k)*13+j), false, trans_rates==0.0)
                                #remove_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(22*13+j), false, trans_rates==0.0)
                            end
                        end
                    end
                    # if i == 300000
                    #     for j in 1:13
                    #         #print(10*13+j, lattice.kinesin[10*13+j])
                    #         attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(24*13+j), false)
                    #         #make_kinesin_like!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(22*13+j))
                    #     end
                    # end
                end
                lengthsreps[rep,:,l,m,n] = mean(eachcol(lengths_protos))
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
#lines!(ax, mean(eachrow(lengthsreps[:,2:length(time),2])), label = "low noise")
for i in 1:1
    for j in 1:1
        # fig = Figure(resolution = (1000, 750))
        # ax = Axis(fig[1,1])
        rate_fixed = rate_list[i]
        #factor_fixed = factors[j]
        lines!(ax, mean(eachcol(vec(sum(clusters[:,:,1,1, i, j], dims = 2)))), label = "$rate_fixed, $j")
#        display(fig)
    end
end        
axislegend(ax, position=:rb)
display(fig)

strng_add = "noise004rate0kinesin100"

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
lbls = ["Two state", "Three state"]
#lines!(ax, mean(eachrow(lengthsreps[:,2:length(time),2])), label = "low noise")
for i in 1:1
    for j in 1:1
        rate_fixed = rate_list[i]
        #factor_fixed = factors[j]
        lines!(ax, collect(0:16/100000:(16*(length(time)-2))/100000), mean(eachrow(lengthsreps[1:reps,2:length(time),1, i, j])), label = lbls[i])
    end
end
#lines!([collect(0:16/1000:(16*(length(time)-2))/1000)[500],collect(0:16/1000:(16*(length(time)-2))/1000)[500]], [113,115], label = "Kinesin added")
#lines!([collect(0:16/1000:(16*(length(time)-2))/1000)[2500],collect(0:16/1000:(16*(length(time)-2))/1000)[2500]], [113,115], label = "Kinesin added")


#lines!(ax, collect(0:16/1000:(16*899)/1000), mean(eachrow(lengthsreps[1:reps,2:length(time),1])), label = "high noise")
axislegend(ax, position=:rt)
ax.xlabel = "Time (ms)"
ax.ylabel = "Length of protofilament section (nm)"
ax.title = "Lengths at different rates"

display(fig)
save(string("plots/length_relaxation_",strng_add,".png"), fig)

CairoMakie.activate!()

for j in 1:1
    fig = Figure(resolution = (1000, 750))
    ax = Axis(fig[1,1])
    plot!(ax, clusters[:,:,1,1,j,1])
    display(fig)
end
save(string("plots/switching_plot_",strng_add,".png"), fig)


for j in 1:1
    fig = Figure(resolution = (1000, 750))
    ax = Axis(fig[1,1])
    plot!(ax, log.(transpose(lattice_rates[:,:,1,1,j,1 ])))
    display(fig)
end
save(string("plots/switching_rates_",strng_add,".png"), fig)



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
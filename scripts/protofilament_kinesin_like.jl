if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


    Nt = 128000
    stp = 8
    time = collect(0:stp:Nt)



    conf = from_toml(MicrotubuleConfig, "config/protofilament_kinesin.toml")

    lattice, bead_info = initialise(conf)




    #attach_kinesin!()


    clusters = falses(length(time), length(lattice.x),1)
    kinesins = falses(length(time), length(lattice.x),1)

    # for i in 1:23
    #     attach_kinesin!(lattice2, bead_info2, conf.lattice.N, conf.lattice.S, Int(i*13 + 4))
    #     attach_kinesin!(lattice2, bead_info2, conf.lattice.N, conf.lattice.S, Int(i*13 + 5))
    # end
    #attach_kinesin!(lattice1, bead_info1, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    #attach_kinesin!(lattice1, bead_info1, conf.lattice.N, conf.lattice.S, Int(5*13 + 5))

    # attach_kinesin!(lattice2, bead_info2, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    # attach_kinesin!(lattice2, bead_info2, conf.lattice.N, conf.lattice.S, Int(5*13 + 5))


    # attach_kinesin!(lattice3, bead_info3, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    # attach_kinesin!(lattice3, bead_info3, conf.lattice.N, conf.lattice.S, Int(5*13 + 5))

    #attach_kinesin!(lattice1, bead_info1, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    #attach_kinesin!(lattice1, bead_info1, conf.lattice.N, conf.lattice.S, Int(5*13 + 5))
   
    #attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(15))

    positions = []



    @showprogress for i in 1:Nt
        iterate!(lattice, bead_info, conf, conf.iter_pars,1)
        if i % stp == 0
            kinesins[Int(i/stp) + 1,:,1] = lattice.kinesin[:,1]
            clusters[Int(i/stp) + 1,:,1] = lattice.kinesin[:,2]
            push!(positions, deepcopy(lattice.x))

        end
    end


CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])

ax.xlabel = "Time stamp"
ax.ylabel = "number of kinesin-bound-like tubulin"
ax.title = "cluster growth, lateral spread"

lines!(ax, vec(sum(clusters[:,:,1], dims = 2)))
display(fig)
save("plots/cluster_growth_energy_proto.png", fig)


fig = Figure(resolution = (1500, 1500))
ax = Axis(fig[1,1], aspect=AxisAspect(100))
limits!(-30,30, -5, 24*5)
frame = Observable(1)
position1 = ([i*4.05 for i in 1:length(lattice.x)])
position2 = ([0 for i in 1:length(lattice.x)])
clrs = @lift([BEADCOLORS[false, clusters[10*$frame, i,1]] for i in 1:length(lattice.x)])

fig = scatter(position1, position2, color = clrs,
    axis = (title = @lift("Frame = $(round(10*$frame, digits = 1))"), aspect=AxisAspect(5), limits = (0,100,-0.5,0.5)), marker = :circle)
#scatter!(position1, position2, color = :red, markersize = 15)
display(fig)

framerate = 60
frames = 1:1600

record(fig, "figures/time_animation_proto.mp4", frames;
        framerate = framerate) do fr
    frame[] = fr
    print(fr)
end


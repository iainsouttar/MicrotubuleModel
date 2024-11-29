if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


    Nt = 12800
    stp = 8
    time = collect(0:stp:Nt)


    conf = from_toml(MicrotubuleConfig, "config/eulerMTkinesinLikestoch.toml")

    lattice, bead_info = initialise(conf)
    lattice1, bead_info1 = initialise(conf)
    lattice2, bead_info2 = initialise(conf)
    lattice3, bead_info3 = initialise(conf)



    print(lattice.rates)

    #attach_kinesin!()


    clusters = falses(length(time), length(lattice.x),4)
    kinesins = falses(length(time), length(lattice.x),4)

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
    attach_kinesin!(lattice1, bead_info1, conf.lattice.N, conf.lattice.S, Int(8*13 + 5))
   
    #attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(5*13 + 4))
    attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(8*13 + 5))
    #attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(6*13 + 5))
    #attach_kinesin!(lattice, bead_info, conf.lattice.N, conf.lattice.S, Int(7*13 + 5))


    positions = []
    positions1 = []
    positions2 = []
    positions3 = []

    noise = 0.004
    @showprogress for i in 1:Nt
        #iterate!(lattice120fourteen2, bead_info120fourteen2, conf120fourteen2, conf120fourteen2.iter_pars)
        iterate!(lattice, bead_info, conf, conf.iter_pars,0, 0.00000000000001, noise)
        #iterate!(lattice1, bead_info1, conf, conf.iter_pars,0,0.0)
        #iterate!(lattice2, bead_info2, conf, conf.iter_pars,0, 0.0)
        #iterate!(lattice3, bead_info3, conf, conf.iter_pars,0, 0.0)


        if i % stp == 0
            kinesins[Int(i/stp) + 1,:,1] = lattice.kinesin[:,1]
            clusters[Int(i/stp) + 1,:,1] = lattice.kinesin[:,2]
            push!(positions, deepcopy(lattice.x))
            kinesins[Int(i/stp) + 1,:,2] = lattice1.kinesin[:,1]
            clusters[Int(i/stp) + 1,:,2] = lattice1.kinesin[:,2]
            push!(positions1, deepcopy(lattice1.x))
            kinesins[Int(i/stp) + 1,:,3] = lattice2.kinesin[:,1]
            clusters[Int(i/stp) + 1,:,3] = lattice2.kinesin[:,2]
            push!(positions2, deepcopy(lattice2.x))
            kinesins[Int(i/stp) + 1,:,4] = lattice3.kinesin[:,1]
            clusters[Int(i/stp) + 1,:,4] = lattice3.kinesin[:,2]
            push!(positions3, deepcopy(lattice3.x))
        end
    end



# Nt = 160000
# stp = 80
# time = collect(0:stp:Nt)



# confstoch = from_toml(MicrotubuleConfig, "config/eulerMTkinesinLikestoch.toml")

# latticestoch, bead_infostoch = initialise(confstoch)
# latticestoch1, bead_infostoch1 = initialise(confstoch)
# latticestoch2, bead_infostoch2 = initialise(confstoch)

# #attach_kinesin!()

# clusters = falses(length(time), length(latticestoch.x),3)
# kinesins = falses(length(time), length(latticestoch.x),3)
# positionsstoch = []
# positionsstoch1 = []
# positionsstoch2 = []

# @showprogress for i in 1:Nt
#     #iterate!(latticestoch120fourteen2, bead_infostoch120fourteen2, conf120fourteen2, conf120fourteen2.iter_pars)
#     iterate!(latticestoch, bead_infostoch, conf, conf.iter_pars,0)
#     iterate!(latticestoch1, bead_infostoch1, conf, conf.iter_pars,1)
#     iterate!(latticestoch2, bead_infostoch2, conf, conf.iter_pars,2)

#     if i % stp == 0
#         kinesins[Int(i/stp) + 1,:,1] = latticestoch.kinesin[:,1]
#         clusters[Int(i/stp) + 1,:,1] = latticestoch.kinesin[:,2]
#         push!(positionsstoch, deepcopy(latticestoch.x))
#         kinesins[Int(i/stp) + 1,:,2] = latticestoch1.kinesin[:,1]
#         clusters[Int(i/stp) + 1,:,2] = latticestoch1.kinesin[:,2]
#         push!(positionsstoch1, deepcopy(latticestoch1.x))
#         kinesins[Int(i/stp) + 1,:,3] = latticestoch2.kinesin[:,1]
#         clusters[Int(i/stp) + 1,:,3] = latticestoch2.kinesin[:,2]
#         push!(positionsstoch2, deepcopy(latticestoch2.x))
#     end
# end


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

f = plot_3d(lattice, bead_info)
f



CairoMakie.activate!()
fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])

ax.xlabel = "Time stamp"
ax.ylabel = "number of kinesin-bound-like tubulin"
ax.title = "cluster growth"

lines!(ax, vec(sum(clusters[:,:,1], dims = 2)), label="higher")
lines!(ax, vec(sum(clusters[:,:,2], dims = 2)), label="lower")
lines!(ax, vec(sum(clusters[:,:,3], dims = 2)), label="higher, no kin")
lines!(ax, vec(sum(clusters[:,:,4], dims = 2)), label="lower, no kin")
axislegend(ax, position=:lc)
display(fig)
save("plots/cluster_growth_energy_kin.png", fig)



CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
plot!(ax, clusters[:,:,1])
display(fig)
save("plots/clusters_plot_energy1.png", fig)


CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
plot!(ax, clusters[:,:,2])
display(fig)
save("plots/clusters_plot_energy2.png", fig)

CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
plot!(ax, clusters[:,:,3])
display(fig)
save("plots/clusters_plot_energy3.png", fig)

CairoMakie.activate!()

fig = Figure(resolution = (1000, 750))
ax = Axis(fig[1,1])
plot!(ax, clusters[:,:,4])
display(fig)
save("plots/clusters_plot_energy4.png", fig)


#create a network
g = SimpleGraph(length(lattice.x))
for i in 1:length(lattice.x)
    for edge in bead_info[i].bonds
        add_edge!(g,i,edge)
    end
end


node_states = Vector(clusters[end,:, 1])
# Compute assortativity based on node states
function clustering_metrics(g::SimpleGraph, node_states::Vector{Bool})
    n = nv(g)
    degrees = [degree(g, i) for i in 1:n]
    total_edges = ne(g)
    sum_edges = 0
    sum_degrees_squared = 0
    num_kin_like = sum(node_states)
    counter = 0

    for edge in edges(g)
        u, v = src(edge), dst(edge)
        sum_edges += node_states[u] * node_states[v]
    end

    metric = sum_edges/num_kin_like
    
    for i in 1:n
        sum_degrees_squared += degrees[i]^2
    end


    for v in vertices(g)[node_states]
        counter += sum(node_states[neighbors(g,v)])/(num_kin_like*length(neighbors(g,v)))
    end
    
    assortativity = (total_edges * sum_edges - sum_degrees_squared) / (total_edges^2 - sum_degrees_squared)
    return assortativity, metric, counter
end

# function clustering_metrics(g::SimpleGraph, node_states::Vector{Bool})
#     n = nv(g)
#     counter = 0
#     counter_edges = 0
#     num_kin_like = sum(node_states)
#     for v in nodes(g)[node_states]
#         counter += sum(node_states[v.neighbours])/(num_kin_like*length(v.neighbours))
#     end

#     for edge in edges(g)
#         u, v = src(edge), dst(edge)
#         counter_edges += node_states[u] * node_states[v]
#     return counter
# end



# assortativity, met, cnter = clustering_metrics(g, node_states)
# println("Assortativity coefficient based on node states: ", assortativity)


# assortativity = zeros(length(clusters[:,1,1]))
# edges_clustering = zeros(length(clusters[:,1,1]))
# prop_edges = zeros(length(clusters[:,1,1]))

# for i in 1:length(clusters[:,1,1])
#     node_states = Vector(clusters[i,:,1])
#     assortativity[i], edges_clustering[i], prop_edges[i] = clustering_metrics(g, node_states)
# end
# CairoMakie.activate!()
# fig = Figure(resolution = (1000, 750))
# ax = Axis(fig[1,1])

# ax.xlabel = "Time stamp"
# ax.ylabel = "metric"
# ax.title = "clusters energy"

# lines!(ax, prop_edges)
# lines!(ax, assortativity)
# lines!(ax, edges_clustering)
# display(fig)
# save("plots/cluster_metrics_energy1.png", fig)

#prop_edges_base = 

# plot(prop_edges.-(vec(sum(clusters[:,:,1], dims = 2)))/312)


# K = 20
# prop_edges_avg = [mean(prop_edges[i:i+K]) for i in 1:(length(prop_edges)-K)]
# assortativity_avg = [mean(assortativity[i:i+K]) for i in 1:(length(prop_edges)-K)]
# edges_clustering_avg = [mean(edges_clustering[i:i+K]) for i in 1:(length(prop_edges)-K)]

# CairoMakie.activate!()
# fig = Figure(resolution = (1000, 750))
# ax = Axis(fig[1,1])

# ax.xlabel = "Time stamp"
# ax.ylabel = "metric"
# ax.title = "clusters"

# lines!(ax, prop_edges_avg)
# lines!(ax, assortativity_avg)
# lines!(ax, edges_clustering_avg)
# display(fig)


# assortativity_avg = [mean(assortativity[i:i+K]) for i in 1:(length(time)-K)]
# plot(assortativity_avg)

ColorSchemes.colorschemes[:nature] = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

COLORSCHEME = ColorScheme([colorant"#E64B35",colorant"#4DBBD5",colorant"#00A087",colorant"#3C5488", colorant"#F39B7F", colorant"#8491B4", colorant"#91D1C2"])

BEADCOLORS = Dict(
    (true, false) => COLORSCHEME.colors[1], 
    (false, false) => COLORSCHEME.colors[2],
    (true, true) => COLORSCHEME.colors[3],
    (false, true) => COLORSCHEME.colors[4]
)

GLMakie.activate!()
GLMakie.closeall()
# Create figure and axes
fig = Figure(resolution = (1000, 1500))
ax = Axis3(fig[1,1])
T=20
# Create scatter plot
for i in 1:length(lattice.x)
    scatter!(ax,positions[T][i], color = BEADCOLORS[false, clusters[T, i,1]], marker = :circle)
    #scatter!(ax, lattice.x, color = :red)
end

# Set title
#fig[1, 1] = title("3D Scatter Plot")
# Show the plot
display(fig)



fig = Figure(resolution = (1500, 1500))
ax = Axis(fig[1,1], aspect=AxisAspect(10))
limits!(-30,30, -5, 24*5)
frame = Observable(1)
per_frame = 1
position1 = @lift([positions[per_frame*$frame][i][3] for i in 1:length(lattice.x)])
position2 = @lift([positions[per_frame*$frame][i][2] for i in 1:length(lattice.x)])
clrs = @lift([BEADCOLORS[false, clusters[per_frame*$frame, i,1]] for i in 1:length(lattice.x)])

fig = scatter(position1, position2, color = clrs,
    axis = (title = @lift("Frame = $(round($frame, digits = 1))"), aspect=AxisAspect(8), limits = (0,110,-30,30)), marker = :circle)
#scatter!(position1, position2, color = :red, markersize = 15)
display(fig)

framerate = 20
frames = 1:Int((length(time)-1)/per_frame)

record(fig, "figures/time_animation_both_spread.mp4", frames;
        framerate = framerate) do fr
    frame[] = fr
    print(fr)
end




fig = Figure(resolution = (1500, 1500))
ax = Axis(fig[1,1], aspect=AxisAspect(10))
limits!(-30,30, -5, 24*5)

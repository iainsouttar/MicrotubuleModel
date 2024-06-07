using Logging
using Parameters: @unpack
using ProgressMeter
using Configurations
using DelimitedFiles


#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

Nt = 100_000
stp = 200
path = "results/raw"

conf = from_toml(MicrotubuleConfig, "config/bending_anim.toml")
conf = set_bond_angles(conf)
lattice, bead_info = MicrotubuleSpringModel.initialise(conf)


filename = "kinesin-$(conf.lattice.num_rings).csv"
filename_bool = "kinesin-bool-$(conf.lattice.num_rings).csv"

data = Matrix{Float64}(zeros(Float64, (length(lattice.x)*3,1)))
for i in 1:length(lattice.x)
    data[3*(i-1)+1:3*i,1] .= lattice.x[i]
end
open(path*"/"*filename, "w") do io
    writedlm(io, data', ',')
end


data = Matrix{Int}(zeros(Int, (length(lattice.x),1)))
for i in 1:length(lattice.x)
    data[i,1] = Int(lattice.kinesin[i])
end
open(path*"/"*filename_bool, "w") do io
    writedlm(io, data', ',')
end

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        MicrotubuleSpringModel.append_to_csv(filename, lattice.x)
        MicrotubuleSpringModel.append_to_csv(filename_bool, Int.(lattice.kinesin))
    end
end

using Accessors

for i in 1:length(lattice.x)
    if i % 13 ∈ (6,7,8,9,10)
        b = bead_info[i]
        lattice.kinesin[i] = true
        pos = get_intra(i, total, b.α)
        bead_info[i] = Accessors.@set b.lin_consts[pos] = 6.0
        bead_info[i] = Accessors.@set b.lengths[pos] = 4.5
    end
end


conf = @set conf.external_force = MicrotubuleSpringModel.NoExternalForce()

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        MicrotubuleSpringModel.append_to_csv(filename, lattice.x)
        MicrotubuleSpringModel.append_to_csv(filename_bool, Int.(lattice.kinesin))
    end
end


for i in 1:length(lattice.x)
    if i % 13 ∈ (6,7,8,9,10)
        b = bead_info[i]
        lattice.kinesin[i] = false
        pos = get_intra(i, total, b.α)
        bead_info[i] = Accessors.@set b.lin_consts[pos] = 3.5
        bead_info[i] = Accessors.@set b.lengths[pos] = 4.05
    end
end

@showprogress for i in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        MicrotubuleSpringModel.append_to_csv(filename, lattice.x)
        MicrotubuleSpringModel.append_to_csv(filename_bool, Int.(lattice.kinesin))
    end
end

################################################################

num_rings = conf.lattice.num_rings
N = num_rings*13

open(path*"/"*filename, "r") do io
    global steps = countlines(io)
end
x = Vector{BeadPos}(undef, N)
kinesin = Vector{Bool}(undef, N)

open(path*"/"*filename_bool) do io_bool
    open(path*"/"*filename) do io
        iterator = CSV.Rows(io, header=false)
        row, state = iterate(iterator)
        update_positions!(x, row, N)
        iterator_bool = CSV.Rows(io_bool, header=false)
        row_bool, _ = iterate(iterator_bool)
        @. kinesin = parse(Int, row_bool)

        x_obs = Observable(x)
        pts = @lift([Point2f(xi[1],xi[3]) for xi in $(x_obs)])

        kin_obs = Observable(kinesin)
        colors = @lift([MicrotubuleSpringModel.COLORS[(iseven(i÷13),k)] for (i,k) in enumerate($(kin_obs))])

        CairoMakie.activate!()
        f = Figure(resolution=(1500,1500), backgroundcolor=colorant"#111111")
        ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
        limits!(-30,50,-10,num_rings*4.05+10)
        hidedecorations!(ax)
        hidespines!(ax)
        scatter!(ax, pts, color=colors, marker=:circle, markersize=40, strokecolor=colorant"#111111")

        p = Progress(steps)
        record(f, "figures/kinesin-$num_rings.mp4", zip(iterator, iterator_bool)) do (row, row_bool)
            update_positions!(x, row, N)
            x_obs[] = x
            @. kinesin = parse(Int, row_bool)
            kin_obs[] = kinesin
            next!(p)
        end
    end
end
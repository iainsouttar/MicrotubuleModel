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

function update_positions!(x::AbstractVector{BeadPos}, row, N)
    for j in 1:N
        x[j] = BeadPos(parse(Float64,row[3*j-2]),parse(Float64,row[3*j-1]), parse(Float64,row[3*j]))
    end
    return
end

function update_kinesin!(x, row, N)
    for j in 1:N
        x[j] = parse(Int,row[j])
    end
    return
end

Nt = 100_000
stp = 200
filename = "bending.csv"
path = "results/raw"

conf = from_toml(MicrotubuleConfig, "config/bending_anim.toml")
conf = set_bond_angles(conf)
beads, bead_info = MicrotubuleSpringModel.initialise(conf)

@unpack num_rings = conf.lattice

data = Matrix{Float64}(zeros(Float64, (length(beads)*3,1)))
for i in 1:length(beads)
    data[3*(i-1)+1:3*i,1] .= beads.x[i]
end
open(path*"/"*filename, "w") do io
    writedlm(io, data', ',')
end

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, conf, conf.iter_pars)
    if i % stp == 0
        MicrotubuleSpringModel.append_to_csv(filename, beads.x)
    end
end

#################################################################

open(path*"/"*filename, "r") do io
    global steps = countlines(io)
end
x = Vector{BeadPos}(undef, N)
kinesin = Vector{Bool}(undef, N)


open(path*"/"*filename) do io
    iterator = CSV.Rows(io, header=false)
    row, _ = iterate(iterator)

    x_obs = Observable(x)
    pts = @lift([Point2f(xi[3],xi[1]) for xi in $(x_obs)])

    CairoMakie.activate!()
    f = Figure(resolution=(1600,900), backgroundcolor=colorant"#111111")
    ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
    limits!(-10,num_rings*4.05+10, -15,num_rings*1.05)
    hidedecorations!(ax)
    hidespines!(ax)
    scatter!(ax, pts, color=colors, marker=:circle, markersize=40, strokecolor=colorant"#111111")
    f

    p = Progress(steps)
    record(f, "figures/bending-$num_rings.mp4", iterator) do row
        update_positions!(x, row, N)
        x_obs[] = x
        next!(p)
    end
end

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using CSV

num_rings = 30
filename = "fluctuations-free-$num_rings.csv"
path = "results/raw"

function update_positions!(x, row, N)
    for j in 1:N
        x[j] = BeadPos(parse(Float64,row[3*j-2]),parse(Float64,row[3*j-1]), parse(Float64,row[3*j]))
    end
    return
end

N = num_rings*13
open(path*"/"*filename, "r") do io
    global steps = countlines(io)
end
x = Vector{BeadPos}(undef, N)

open(path*"/"*filename) do io
    iterator = CSV.Rows(io, header=false)
    row, state = iterate(iterator)
    update_positions!(x, row, N)

    x_obs = Observable(x)
    pts = @lift([Point2f(xi[1],xi[3]) for xi in $(x_obs)])

    CairoMakie.activate!()
    f = Figure(resolution=(1500,1500), backgroundcolor=colorant"#F8F5EE")
    ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#F8F5EE")
    limits!(-15,15,-10,num_rings*4.05+10)
    hidedecorations!(ax)
    hidespines!(ax)
    scatter!(ax, pts, color=[MicrotubuleSpringModel.COLORS[(iseven(i÷13),false)] for i in 1:N], marker=:circle, markersize=50, strokecolor=:black)
    f

    p = Progress(steps)
    record(f, "figures/flat-$num_rings-light.mp4", iterator) do row
        update_positions!(x, row, N)
        x_obs[] = x
        next!(p)
    end
end
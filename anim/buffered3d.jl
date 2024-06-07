"""
Animates a 3d visualisation of the MT given raw data from a csv file.
Loads the data a row at a time in case the dataset is large.
"""

#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


filename = "fluctuations-free-30.csv"
path = "results/raw"

function update_positions!(x, row, N)
    for j in 1:N
        x[j] = BeadPos(parse(Float64,row[3*j-2]),parse(Float64,row[3*j-1]), parse(Float64,row[3*j]))
    end
    return
end

N = 10
open(path*"/"*filename, "r") do io
    global steps = countlines(io)
end
print(steps)
x = Vector{BeadPos}(undef, N)

print(N)
println(N)

open(path*"/"*filename) do io
    iterator = CSV.Rows(io, header=false)
    row, state = iterate(iterator)
    update_positions!(x, row, N)

    GLMakie.activate!()
    GLMakie.closeall()
    scene = Scene(resolution=(1200,900), backgroundcolor=colorant"#111111")
    cam3d!(scene)
    bead_plots = plot_individual!(scene, x, a=8.0)
    GLMakie.scale!(scene, 0.03, 0.03, 0.03)
    center!(scene)

    p = Progress(steps)
    Makie.record(scene, "figures/anim-free-100.mp4", iterator) do row
        update_positions!(x, row, N)
        for (n, x_i) in enumerate(x)
            translate!(bead_plots[n], x_i...)
        end
        GLMakie.scale!(scene, 0.05, 0.05, 0.05)
        center!(scene)
        next!(p)
    end
end
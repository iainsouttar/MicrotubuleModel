using Logging
using Parameters: @unpack
using ProgressMeter
using Configurations


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
filename = "fluctuations.csv"
path = "results/raw"

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/stochastic.toml")
conf = set_bond_angles(conf)
beads, bead_info = MicrotubuleSpringModel.initialise(conf)

data = Matrix{Float64}(zeros(Float64, (length(beads.x)*3,1)))
for i in 1:lastindex(beads.x)
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

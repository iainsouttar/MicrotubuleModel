if isdefined(@__MODULE__, :LanguageServer)
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


conf = from_toml(PatchConfig, "config/patch.toml")
conf = MicrotubuleSpringModel.set_bond_angles(conf)

lattice, dirs = MicrotubuleSpringModel.initialise(conf)

# @showprogress for i in 1:10000
#     iterate!(lattice, conf, dirs)
# end

GLMakie.closeall()
scene = plot(lattice, dirs)
scene
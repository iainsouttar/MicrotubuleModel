if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end


function main!(conf, Nt)

    beads, bead_info, dirs = burnin(conf, 1000)

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    end

    for (i,b) in enumerate(beads)
        if i % 13 ∈ (4,5,6,7,8,9)
            b.kinesin = true
        end
    end

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    end

    conf = @set conf.external_force = MicrotubuleSpringModel.NoExternalForce()

    @showprogress for i in 1:Nt
        iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    end

    return beads, bead_info
end


Nt = 10_000


conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/bending_stiffness.toml")
conf = set_bond_angles(conf)

beads, bead_info, dirs = burnin(conf, 1000)

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
end

for (i,b) in enumerate(beads)
    if i % 13 ∈ (4,5,6,7,8,9)
        b.kinesin = true
    end
end

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
end

conf = @set conf.external_force = MicrotubuleSpringModel.NoExternalForce()

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
end


CairoMakie.activate!()
f = Figure(resolution=(1000,1500))
ax = Axis(f[1,1], aspect=DataAspect())
plot_flat!(ax, beads, bead_info)
hidedecorations!(ax)
hidespines!(ax)
f


GLMakie.activate!()
GLMakie.closeall()
scene = plot(beads, bead_info)
scene

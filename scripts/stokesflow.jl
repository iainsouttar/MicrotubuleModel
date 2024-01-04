#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

#####################################################

conf = from_toml(MicrotubuleSpringModel.RotationConfig, "config/youngs_modulus.toml")

conf = set_bond_angles(conf)

beads, bead_info, dirs = MicrotubuleSpringModel.initialise(conf)

###################################################################
beads = deepcopy(beads_cpy)
bead_info = deepcopy(bead_info_cpy)


conf = @set conf.external_force = ym_conf


Nt = 50_000
step = 50
time = 0:step:Nt
ext = zeros(Nt÷step+1)

L0 = microtubule_length(beads, conf.lattice)

@showprogress for i in 1:Nt
    iterate!(beads, bead_info, dirs, conf, conf.iter_pars)
    if i % step == 0
        ext[i÷step+1] = microtubule_length(beads, conf.lattice) - L0
    end
end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, time, ext)
f


function microtubule_length(beads, consts)
    @unpack N, num_rings = consts
    Ntot = lastindex(beads)
    d = 0.0
    for i in 1:N
        d += norm(beads[Ntot-N+i].x - beads[i].x)
    end
    return d/N
end
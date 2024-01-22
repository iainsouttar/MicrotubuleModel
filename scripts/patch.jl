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

Nt = 50_000
step = 200
path = "results/raw"

conf = from_toml(MicrotubuleConfig, "config/kinesin.toml")
conf = set_bond_angles(conf)
lattice, bead_info = MicrotubuleSpringModel.initialise(conf)

l0_kin = conf.spring_consts.l0_in_kin
k_kin = conf.spring_consts.k_in_kin
num_rings = conf.lattice.num_rings
N = num_rings*13

using Accessors

i = 30
lattice.kinesin[i] = true
b = bead_info[i]
pos = get_intra(i, N, b.α)
bead_info[i] = Accessors.@set b.lin_consts[pos] = k_kin
bead_info[i] = Accessors.@set b.lengths[pos] = l0_kin
j = bead_info[i].bonds[pos]
b = bead_info[j]
lattice.kinesin[j] = true
pos = get_intra(j, total, b.α)
bead_info[j] = Accessors.@set b.lin_consts[pos] = k_kin
bead_info[j] = Accessors.@set b.lengths[pos] = l0_kin


n_bonds = sum(length(b.bonds)*3 for b in bead_info)
pts = Vector{Point2f}(undef, n_bonds)
ext = Vector{Float64}(undef, n_bonds)

function extensions!(pts, ext, positions, info)
    idx = 1
    for (x1, b) in zip(positions, info)
        for (j,(y, l)) in enumerate(zip(b.bonds, b.lengths))
            x2 = positions[y]
            ext[idx:idx+2] .= (norm(x2 .- x1) - l)/l
            pts[idx] = Point2f(x1[1],x1[3])
            pts[idx+1] = Point2f(x2[1],x2[3])
            pts[idx+2] = Point2f(NaN,NaN)
            idx += 3
        end
    end
    return pts, ext
end


for t in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
end

latt_obs = Observable(lattice)
xs = @lift(Vector{Point2f}([Point2f(xi[1],xi[3]) for xi in $(latt_obs).x]))
colors = @lift([MicrotubuleSpringModel.COLORS[(iseven(i÷13),k)] for (i,k) in enumerate($(latt_obs).kinesin)])

extensions!(pts, ext, lattice.x, bead_info)
pts_obs, ext_obs = Observable(pts), Observable(ext)
colorrange = @lift(extrema($(ext_obs)))
#colorrange = (-0.05,0.05)

CMAP = :plasma
f = Figure(resolution=(1200,1000), backgroundcolor=colorant"#111111")
ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
limits!(-12,12,-2,num_rings*4.05+10)
hidedecorations!(ax)
hidespines!(ax)#colorrange=colorrange)
# scatter!(
#     ax, xs, 
#     color=colors, marker=:circle, markersize=100, strokecolor=colorant"#111111", strokewidth=3
# )
lines!(ax, pts_obs, color=ext_obs, linewidth=15, colormap=Reverse(CMAP), colorrange=colorrange)

Colorbar(
    f[1,2], 
    label=L"\text{Strain } (\Delta L/L_0)", labelcolor=:white, limits = colorrange, width=50, height=900, colormap = Reverse(CMAP), ticklabelcolor=:white, ticklabelsvisible=true
)
f

save("figures/kinesin-strain-3.png", f, px_per_unit=2)

p = Progress(Nt÷step+1)
frames = range(1, Nt÷step+1, step=1)
record(f, "figures/single-kinesin.mp4", frames, framerate=25) do i
    for t in 1:step
        iterate!(lattice, bead_info, conf, conf.iter_pars)
    end
    extensions!(pts, ext, lattice.x, bead_info)
    latt_obs[] = lattice
    pts_obs[] = pts
    ext_obs[] = ext
    next!(p)
end



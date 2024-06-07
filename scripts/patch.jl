"""
This script creates a short MT then changes the bond strength and length for a single dimer. It then visualises the resulting strains of the lattice once the new equilibrium has been found.

"""

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

using Accessors

"""Change the constants for a dimer, ie. an alpha+beta"""
function attach_kinesin!(lattice, bead_info, N, S, i::Int)
    #we want this to change a dimer, rather than between dimers, this
    #is what the pos business is doing
    lattice.kinesin[i,:] = [true, true]
    b = bead_info[i]
    pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)


    #bead_info[i] = Accessors.@set b.lin_consts[pos] = k_kin
    #bead_info[i] = Accessors.@set b.lin_consts[pos] = k_kin
    #bead_info[i] = Accessors.@set b.bend_consts[pos] = K_in_kin
    #bead_info[i] = Accessors.@set b.bend_consts[pos] = 0.0

    bead_info[i] = Accessors.@set b.lengths[pos] = l0_kin
    #bead_info[i] = Accessors.@set b.directions[1] = [0.0,0.0,1.0]

    j = bead_info[i].bonds[pos]
    b = bead_info[j]
    lattice.kinesin[j,:] = [true, true]
    pos = intra_dimer_index(j, length(lattice.x), b.α, N, S)
    #print(j,pos)



    #bead_info[j] = Accessors.@set b.lin_consts[pos] = k_kin
    bead_info[j] = Accessors.@set b.lengths[pos] = l0_kin
    #bead_info[j] = Accessors.@set b.bend_consts[pos] = K_in_kin

    #bead_info[j] = Accessors.@set b.bend_consts[pos] = 0.0
end


"""Change the constants for a dimer, ie. an alpha+beta"""
function make_kinesin_like!(lattice, bead_info, N, S, i::Int)
    #we want this to change a dimer, rather than between dimers, this
    #is what the pos business is doing
    lattice.kinesin[i,2] = true
    b = bead_info[i]
    pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)


    bead_info[i] = Accessors.@set b.lin_consts[pos] = k_kin
    #bead_info[i] = Accessors.@set b.lin_consts[pos] = k_kin
    bead_info[i] = Accessors.@set b.bend_consts[pos] = K_in_kin
    bead_info[i] = Accessors.@set b.bend_consts[pos] = 0.0

    bead_info[i] = Accessors.@set b.lengths[pos] = l0_kin
    #bead_info[i] = Accessors.@set b.directions[1] = [0.0,0.0,1.0]



    bead_info[j] = Accessors.@set b.lin_consts[pos] = k_kin
    bead_info[j] = Accessors.@set b.lengths[pos] = l0_kin
    bead_info[j] = Accessors.@set b.bend_consts[pos] = K_in_kin

    bead_info[j] = Accessors.@set b.bend_consts[pos] = 0.0
end

"""Calculate extensions ready for plotting bonds"""
function extensions!(pts, ext, positions, info)
    idx = 1
    for (x1, b) in zip(positions, info)
        for (j,(y, l)) in enumerate(zip(b.bonds, b.lengths))
            x2 = positions[y]
            ext[idx:idx+2] .= (norm(x2 .- x1) - l)/l
            pts[idx] = Point2f(x1[1],x1[3])
            pts[idx+1] = Point2f(x2[1],x2[3])
            # put nans between bonds for discontinuous lines
            pts[idx+2] = Point2f(NaN,NaN)
            idx += 3
        end
    end
    return pts, ext
end

Nt = 50_00
stp = 200
path = "results/raw"

conf = from_toml(MicrotubuleConfig, "config/eulerMTkinesinLike.toml")
conf = set_bond_angles(conf)
lattice, bead_info = MicrotubuleSpringModel.initialise(conf)

l0_kin = conf.spring_consts.l0_in_kin
k_kin = conf.spring_consts.k_in_kin
K_in_kin = conf.spring_consts.K_in_kin

attach_kinesin!(lattice, bead_info,conf.lattice.N, conf.lattice.S, 13*12+4)
#attach_kinesin!(lattice, bead_info,conf.lattice.N, conf.lattice.S, 13*12+5)
#attach_kinesin!(lattice, bead_info,conf.lattice.N, conf.lattice.S, 13*12+3)
attach_kinesin!(lattice, bead_info,conf.lattice.N, conf.lattice.S, 13*14+4)


@showprogress for t in 1:Nt
    iterate!(lattice, bead_info, conf, conf.iter_pars)
end

#####################################################

# Visualise stretched bonds
num_rings = 24

n_bonds = sum(length(b.bonds)*3 for b in bead_info)
pts = Vector{Point2f}(undef, n_bonds)
ext = Vector{Float64}(undef, n_bonds)

xs = Vector{Point2f}([Point2f(xi[1],xi[3]) for xi in lattice.x])
#colors = [MicrotubuleSpringModel.COLORS[(iseven(i÷13),k)] for (i,k) in enumerate(lattice.kinesin)]

extensions!(pts, ext, lattice.x, bead_info)

CMAP = :plasma
f = Figure(resolution=(1200,1000), backgroundcolor=colorant"#111111")
ax = Axis(f[1,1], aspect=DataAspect(), backgroundcolor=colorant"#111111")
limits!(-12,12,35,18*4.05)
hidedecorations!(ax)
hidespines!(ax)
scatter!(
    ax, xs, 
    marker=:circle, markersize=20, strokecolor=colorant"#111111", strokewidth=3
) #color =colors
lines!(ax, pts, color=ext, linewidth=15, colormap=Reverse(CMAP))

Colorbar(
    f[1,2], 
    label=L"\text{Strain } (\Delta L/L_0)", 
    labelcolor=:white, width=50, height=900,
    colormap = Reverse(CMAP), limits=extrema(ext),
    ticklabelcolor=:white, ticklabelsvisible=true
)
f
display(f)

save("figures/kinesin-strain1.png", f, px_per_unit=2)




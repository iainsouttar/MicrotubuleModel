#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

@testset "Testing lattice creation" begin
    dirs = Dict(
        true => MicrotubuleSpringModel.bond_directions(MicrotubulesSpringModel.AlphaConfirm()),
        false => MicrotubuleSpringModel.bond_directions(MicrotubulesSpringModel.BetaConfirm())
    )
    a = 4.05
    δx = 5.13
    lattice, bead_info = create_lattice(dirs, MicrotubuleSpringModel.SpringConst(), 5, a, δx)
    Ntot = length(lattice)

    @test length(lattice) == 5*13
    @test bead_info[5].α == false
    @test bead_info[23].α == true
    @test bead_info[3*13+5].α == true

    # test that the bonded beads are contained in the closest four beads in the lattice
    pos = [b.x for b in lattice]
    distances = reshape([dot(x1-x2,x1-x2) for x1 in pos for x2 in pos],(Ntot,Ntot))
    distances[distances .== 0] .= Inf
    setdiffs = 0
    for i in 1:Ntot
        lat, (intra,long) = neighbours(i, Ntot)
        neigh = Set((lat...,intra,long))
        delete!(neigh, 0)
        nn = Set(sortperm(distances[:,i])[1:4])
        setdiffs += length(setdiff(neigh,nn))
    end
    @test setdiffs == 0
end
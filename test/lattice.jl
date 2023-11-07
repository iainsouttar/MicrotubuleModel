

@testset "Testing lattice creation" begin
    a = 4.05
    δx = 5.13
    lattice = create_lattice(5, a, δx)
    Ntot = length(lattice)

    @test length(lattice) == 5*13
    @test lattice[5].α == true
    @test lattice[23].α == false
    @test lattice[3*13+5].α == false
    
    @test lattice[14].θ == BeadAngle(0,0,0)

    pos = [b.x for b in lattice]
    distances = reshape([dot(x1-x2,x1-x2) for x1 in pos for x2 in pos],(Ntot,Ntot))
    distances[distances .== 0] .= Inf
    setdiffs = 0
    for i in 1:Ntot
        lat, long = neighbours(i, Ntot)
        neigh = Set((lat...,long...))
        nn = Set(sortperm(distances[:,i])[1:4])
        setdiffs += length(setdiff(neigh,nn))
    end
    @test setdiffs == 0
end


@testset "Testing orientation" begin
    v = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π/4,π/4))
    @test norm(v) ≈ 1.0

    # test confirmation

    a = 4.05
    δx = 5.13
    lattice, bead_info = create_lattice(5, a, δx)

    v, θ = MicrotubuleSpringModel.to_axis_angle(lattice[1].q)

end
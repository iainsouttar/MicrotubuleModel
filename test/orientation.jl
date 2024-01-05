

@testset "Testing orientation" begin
    v = MicrotubuleSpringModel.direc_from_angles(BondAngle(π/4,π/4))
    @test norm(v) ≈ 1.0

    #v_ = transform_orientation()

    theta, thetahat, nhat = MicrotubuleSpringModel.bond_angle(v_prime, normalize(BondDirec(0.5,0.2,1)))

    @test norm(thetahat) ≈ 1.0
    @test norm(nhat) ≈ 1.0

    # test confirmation

    v, θ = MicrotubuleSpringModel.to_axis_angle(lattice[1].q)

    F = zeros(Float64, (3,length(lattice)))
    torque = MicrotubuleSpringModel.angular_forces!(F, lattice[1], lattice, dirs[lattice[1].α], conf.spring_consts.K)
end

@option "alpha struct" struct AlphaConfirm
    north::SVector{2,Float64} = [π/2, -0.2]
    east::SVector{2,Float64} = [π/13, π/2-0.1819]
    south::SVector{2,Float64} = [π/2, π]
    west::SVector{2,Float64} = [π-π/13,π/2+0.1819]
end

@option "beta struct" struct BetaConfirm
    north::SVector{2,Float64} = [π/2, 0]
    east::SVector{2,Float64} = [π/13, π/2-0.1819]
    south::SVector{2,Float64} = [π/2, π+0.2]
    west::SVector{2,Float64} = [π-π/13,π/2+0.1819]
end


function calc_natural_angles(S, N, dx, a)
    r = S*a/N
    R = N*dx/2π
    l = 2*R*sin(π/N)
    ϕ = atan(r,l)
    #δ = 0.2
    δ = 0.0

    α = MicrotubuleSpringModel.AlphaConfirm(
        [π/2, -δ],
        [π/13, π/2-ϕ],
        [π/2, π],
        [π-π/13,π/2+ϕ]
    )
    β = MicrotubuleSpringModel.BetaConfirm(
        [π/2, 0],
        [π/13, π/2-ϕ],
        [π/2, π+δ],
        [π-π/13,π/2+ϕ]
    )
    return α, β
end

"""
`AlphaConfirm` specifies the bond directions for each of the four bonds for an alpha tubulin using an Euler angle representation.
"""
@option "alpha struct" struct AlphaConfirm
    north::SVector{2,Float64} = [-π/2, 0.2]
    east::SVector{2,Float64} = [π/13, π/2+0.1819]
    south::SVector{2,Float64} = [-π/2, π]
    west::SVector{2,Float64} = [π-π/13,π/2-0.1819]
end

"""
`BetaConfirm` specifies the bond directions for each of the four bonds for an beta tubulin using an Euler angle representation.
"""
@option "beta struct" struct BetaConfirm
    north::SVector{2,Float64} = [-π/2, 0]
    east::SVector{2,Float64} = [π/13, π/2+0.1819]
    south::SVector{2,Float64} = [-π/2, π-0.2]
    west::SVector{2,Float64} = [π-π/13,π/2-0.1819]
end


"""
    bond_directions(thetas::Union{AlphaConfirm,BetaConfirm})

Converts bond directions from Euler angle representations to a Vector of Vectors, one for each bond.
"""
function bond_directions(thetas::Union{AlphaConfirm,BetaConfirm})
    vs = [direc_from_angles(BondAngle(t)) for t in [thetas.north,thetas.east,thetas.south,thetas.west]]
    return vs
end

"""
    calc_natural_angles(S, N, dx, a)

Initialise bond directions based on the lattice structure such that the lateral bonds directions match their current state.

# Arguments
- `S::Int`: helix-start number (longitudinal increase between first and last in a ring)
- `N::Int`: protofilament number (number of beads in each ring)
- `dx::Real`: axial distance between adjacent beads in a ring
- `a::Real`: longitudinal distance between adjacent beads in a protofilament
"""
function calc_natural_angles(S, N, dx, a)
    r = S*a/N  # subunit rise
    R = N*dx/2π  # radius of cylinder
    l = 2*R*sin(π/N)  # distance between adjacent beads
    ϕ = atan(r,l)
    δ = 0.2

    α = MicrotubuleSpringModel.AlphaConfirm(
        [π/2, -δ],
        [π/13, π/2+ϕ],
        [π/2, π],
        [π-π/13,π/2-ϕ]
    )
    β = MicrotubuleSpringModel.BetaConfirm(
        [π/2, 0],
        [π/13, π/2+ϕ],
        [π/2, π+δ],
        [π-π/13,π/2-ϕ]
    )
    return α, β
end
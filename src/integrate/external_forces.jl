"""
Add external forces applied to specific beads in the lattice for various experiments.
Default is no external force.

"""


function external_forces!(F, beads, bead_info, consts::IsotropicForce)
    force = SVector{3,Float64}(consts.F, 0, 0)
    @inbounds @fastmath for i in 1:lastindex(beads)
        F[:,i] .+= force
    end
end

function external_forces!(F, beads, bead_info, consts::YoungsModulusTest)
    force = SVector{3,Float64}(0, 0, consts.F)
    Ntot = lastindex(beads)
    @inbounds @fastmath for i in 1:Ntot
        F[:,i] .+= force*(i>Ntot-consts.N-1)
    end
end

function external_forces!(F, beads, bead_info, consts::BendingStiffnessTest)
    force = SVector{3,Float64}(consts.F, 0, 0)
    Ntot = lastindex(beads)
    @inbounds @fastmath for i in 1:Ntot
        F[:,i] .+= force*(i>Ntot-consts.N-1)
    end
end

function external_forces!(F, beads, bead_info, consts::NoExternalForce) end
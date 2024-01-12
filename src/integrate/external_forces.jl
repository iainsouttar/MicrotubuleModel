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
    Ntot = length(beads)
    @inbounds @fastmath for i in Ntot-consts.N-1:Ntot
        F[3,i] += consts.F
    end
end

function external_forces!(F, beads, bead_info, consts::YoungsModulusPF)
    F[3,end] += consts.F
end

function external_forces!(F, beads, bead_info, consts::BendingStiffnessTest)
    Ntot = length(beads)
    @inbounds @fastmath for i in Ntot-consts.N-1:Ntot
        F[1,i] += consts.F
    end
end

function external_forces!(F, beads, bead_info, consts::NoExternalForce) end
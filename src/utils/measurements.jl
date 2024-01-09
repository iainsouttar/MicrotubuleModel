

function microtubule_length(beads, consts)
    N = consts.N
    Ntot = lastindex(beads)
    last = @view beads[Ntot-N:end]
    # average axial distance between first and last rings of beads
    tot = mapreduce((a,b)->abs(a.x[3]-b.x[3]), +, beads[1:N], last)
    return tot / N
end


function deflection_end(beads, original)
    N = length(original)
    Ntot = lastindex(beads)
    last = @view beads[Ntot-N:end]
    # average transverse distance between the end ring and its original position
    tot = mapreduce((a,b)->abs(a.x[1]-b.x[1]), +, original, last)
    return tot / N
end


surface_area(R, t) = 2π*R*t

"""
returns youngs modulus measured in GPa if the Force is nN while the area is nm^2 and the lengths are in the same units.
"""
function youngs_modulus(F, area, extension, L0)
    stress = F / area
    strain = extension / L0
    return stress / strain
end


stiffness(F, L0, f) = F*L0^3/(3*f)

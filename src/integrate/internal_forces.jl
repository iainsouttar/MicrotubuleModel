
using Accessors
"""
    internal_forces_and_torques!(F, F_, torque, lattice, bead_info, consts)

Calculate then aggregate the forces and torques on each bead due to internal interactions.
"""
@inline function internal_forces_and_torques!(F, F_, torque, lattice, bead_info, flag, N, S,K_r,K_r1, Tk_B)
    Ntot = length(lattice.x)

    @inbounds @fastmath @threads for i in 1:Ntot
        b = bead_info[i]
        bonds = lattice.x[b.bonds]
        bond_orients = lattice.q[b.bonds]
        torque[:,i], F[:,i], F_[:,:,i] = bending_and_spring_forces(lattice.x[i], lattice.q[i], bonds, b, bond_orients, i)
        rate = lattice.kinesin[i,2] ? K_r : K_r1

        #flags: 0 is based on kramers reaction theory
        #flags: 1 is based on longitudinal spread of clusters
        #flags: 2 is a constant switching rate
        #flags: 3 is based on lateral spread of clusters
        if flag == 0
            if rate == 0 lattice.rates[i] = 0
            else lattice.rates[i] = rate*exp(calc_switch_rate(N, lattice.x[i], lattice.q[i], bonds, b, bond_orients, i, N, S)/Tk_B) #+ 0.00001*lattice.kinesin[i,2] + 0.0*!lattice.kinesin[i,2]
            end
        elseif flag == 1
            indices = bond_indices(i, Ntot, b.α, N, S)
            lattice.rates[i] = 0.01*sum(lattice.kinesin[b.bonds[indices.!=1], 2] .!= lattice.kinesin[i,2]) + 0.0004*lattice.kinesin[i,2] + 0.0001*!lattice.kinesin[i,2]
        elseif flag == 2
            lattice.rates[i] = 0.1
        elseif flag == 3
            indices = bond_indices(i, N, b.α, N, S)
            lattice.rates[i] = 0.005*sum(lattice.kinesin[b.bonds[indices.==1], 2] .!= lattice.kinesin[i,2])


        end
    end
    @threads for i in 1:Ntot
        b = bead_info[i]
        for j in 1:lastindex(b.bonds)
            @. F[:, b.bonds[j]] += F_[:, j, i]
        end
    end
end



"""Change the constants for a dimer, ie. an alpha+beta"""
#This needs a bit of a tidy up and should allow for different
#conformational state changes. In particular there should be at
#least 4 options, corresponding to the proposed conformational
#changes described in the report.
function make_kinesin_like!(lattice, bead_info, N, S, i::Int)
    #we want this to change a dimer, rather than between dimers, this
    #is what the pos stuff is doing
    if lattice.kinesin[i,2]
        if ~(lattice.kinesin[i,1])
            lattice.kinesin[i,2] = false
            b = bead_info[i]
            pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)

            b = bead_info[i]
            pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)
            lat_bonds = MicrotubuleSpringModel.lat_indices(i, length(lattice.x),b.α, N, S)
            bnds = MicrotubuleSpringModel.bond_indices(i, length(lattice.x),b.α, N, S)
            idxs = findall(x -> x==1, bnds)
            if lat_bonds[1] != 0
                #b = Accessors.@set b.bend_consts[idxs[1]] = 1.0
                #print(bead_info[i].bend_consts)
    
                bead_info[i] = Accessors.@set b.directions[idxs[1]] = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π/13, π/2+0.1819))
                                
                #bead_info[i] = Accessors.@set b.lengths[idxs[1]] = 5.21
            end
            if lat_bonds[2] != 0
                c = bead_info[i]
                c = Accessors.@set c.bend_consts[1+(lat_bonds[1] != 0)] = 1.0
                #print(bead_info[i].bend_consts)
    
                bead_info[i] = Accessors.@set c.directions[idxs[1+(lat_bonds[1] != 0)]] = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π-π/13, π/2-0.1819))
                

                #bead_info[i] = Accessors.@set b.directions[idxs[1+(lat_bonds[1] != 0)]] = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π-π/13, π/2-0.1819))
                #bead_info[i] = Accessors.@set b.lengths[idxs[1+(lat_bonds[1] != 0)]] = 5.21
            end


            #bead_info[i] = Accessors.@set b.lengths[pos] = 4.05
            #bead_info[i] = Accessors.@set b.directions[1] = [0.0,0.0,1.0]

            # j = bead_info[i].bonds[pos]
            # b = bead_info[j]
            # lattice.kinesin[j,2] = false
            # pos = intra_dimer_index(j, length(lattice.x), b.α, N, S)
            #print(j,pos)



            #bead_info[j] = Accessors.@set b.lengths[pos] = 4.05
        end

        # if lattice.kinesin[i,1]
        #     lattice.kinesin[i,1] = false
        #     lattice.kinesin[i,2] = true

        #     b = bead_info[i]
        #     pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)
    
    
        #     bead_info[i] = Accessors.@set b.lengths[pos] = 4.05
        #     #bead_info[i] = Accessors.@set b.directions[1] = [0.0,0.0,1.0]
    
        #     j = bead_info[i].bonds[pos]
        #     b = bead_info[j]
        #     lattice.kinesin[j,1] = false
        #     lattice.kinesin[j,2] = true

        #     pos = intra_dimer_index(j, length(lattice.x), b.α, N, S)
        #     #print(j,pos)
    
    
    
        #     bead_info[j] = Accessors.@set b.lengths[pos] = 4.05
        # end

    else
        lattice.kinesin[i,2] = true
        b = bead_info[i]
        pos = intra_dimer_index(i, length(lattice.x), b.α, N, S)
        lat_bonds = MicrotubuleSpringModel.lat_indices(i, length(lattice.x),b.α, N, S)
        bnds = MicrotubuleSpringModel.bond_indices(i, length(lattice.x),b.α, N, S)
        idxs = findall(x -> x==1, bnds)
        if lat_bonds[1] != 0
            #b = Accessors.@set b.bend_consts[idxs[1]] = 10.0
            #print(bead_info[i].bend_consts)

            bead_info[i] = Accessors.@set b.directions[idxs[1]] = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π/3, π/2+0.1819))
            

            #bead_info[i] = Accessors.@set b.lengths[idxs[1]] = 5.21
        end
        if lat_bonds[2] != 0
            c = bead_info[i]
            #c = Accessors.@set c.bend_consts[idxs[1+(lat_bonds[1] != 0)]] = 10.0
            #print(bead_info[i].bend_consts)

            bead_info[i] = Accessors.@set c.directions[idxs[1+(lat_bonds[1] != 0)]] = MicrotubuleSpringModel.direc_from_angles(MicrotubuleSpringModel.BondAngle(π-π/3, π/2-0.1819))
            #bead_info[i] = Accessors.@set b.lengths[idxs[1+(lat_bonds[1] != 0)]] = 5.21
        end


        ####Now only affects the tubulin it is attached to
        # bead_info[i] = Accessors.@set b.lengths[pos] = 4.25
        # #bead_info[i] = Accessors.@set b.directions[1] = [0.0,0.0,1.0]

        # j = bead_info[i].bonds[pos]
        # b = bead_info[j]
        # lattice.kinesin[j,2] = true
        # pos = intra_dimer_index(j, length(lattice.x), b.α, N, S)
        # #print(j,pos)



        # bead_info[j] = Accessors.@set b.lengths[pos] = 4.25
    end

end



"""
    calc_switch_rate()
"""
function calc_switch_rate(Ntot,x,q,bonds,b,bond_orients,i, N, S)
    indices = bond_indices(i, Ntot, b.α, N, S)
    engy = bead_energy(x, q, indices, bonds, b, bond_orients)
    lat_engy = engy[1] + engy[4]
    long_engy = engy[2] + engy[3] + engy[5] + engy[6]
    energy = lat_engy + long_engy 
    return energy
end


"""
    bending_and_spring_forces(x, q, bonds, b, K)

Calculate 3D torque and force acting on bead `b1` and its neighbours due to the bond angle bending and spring force at `b1`. Updates overall force vectors and returns torque on `b1`.
"""
function bending_and_spring_forces(x, q, bonds, b, bond_orients, idx)
    @unpack lin_consts, bend_consts, lengths, directions, torque_consts, pos_consts = b
    F_ = MMatrix{3, 4, Float64}(undef)
    torque1 = MVector{3,Float64}(0,0,0)
    F1 = MVector{3,Float64}(0,0,0)
    e2 = SVector{3,Float64}(1,0,0)
    #print(x, bonds, "\n")

    #Force constant for being close to the x axis
    K_pos = pos_consts


    for (i, (k, K, l0, dir, bx, bq, K_torque)) in enumerate(zip(lin_consts, bend_consts, lengths, directions, bonds, bond_orients, torque_consts))
        rhat, d = norm_and_mag(bx-x)
        # transform bond direction according to bead orientation
        v = orientate_vector(dir, sign(q))
        # torque from diff between rest direction v and actual r
        τ, F = bending_torque_and_force(v, rhat, d, K)

        if(K_torque ==0)
            τ1 = zeros(3)
            #print("no torsion")
        else
            #torque from torsion
            #normal_vec = normalize(cross(rhat, e2))
            orientationi = orientate_vector(e2, sign(q))
            orientationj = orientate_vector(e2, sign(bq))
            # orientationiproj = orientationi - dot(orientationi, rhat)*rhat
            # orientationjproj = orientationj - dot(orientationj, rhat)*rhat
            # orientationiproj = normalize(orientationiproj)
            # orientationjproj = normalize(orientationjproj)

            #orientationi = orientate_vector(BondDirec(1.0,0.0,0.0), sign(q))
            #orientationj = orientate_vector(BondDirec(1.0,0.0,0.0), sign(bq))
            τ1 = torsion_torque(rhat, d, orientationi, orientationj, K_torque)
            # if (norm(τ1) > 2)
            #     print("\n", orientationi, orientationiproj, orientationj, orientationjproj, "\n")
            # end
        end


        @. F_[:,i] = -F
        torque1 += τ + τ1
        F1 += F
        F1 += spring_force(rhat, d, l0, k)
    end
    #the following might be able to be done by just calculating the norm of the first two coordinates of
    #the position vector x
    #D = cross(x.-(0.0,0.0,1.0), MVector{3,Float64}(0.0,0,1.0))
    if b.α
        D = K_pos*MVector{3,Float64}(-x[1],-x[2],0.0)
        F1 += D
    end
    return torque1, F1, F_
end


"""
    torque_and_force(v, rhat, d, K)

Calculate the torque and force due to the bond angle bending between the natural direction `v` and the actual direction `rhat` with length `d`.
"""
@inline @fastmath function bending_torque_and_force(v, rhat, d, K)


    rdotv = clamp(dot(rhat,v), -1, 1)

    if ((1-rdotv) < 1e-5)
        return zeros(3), zeros(3)
    end

    # τ = K*sin(θ)*n̂
    torque = K*sqrt(1.0-rdotv^2)*normalize(cross(rhat, v))
    # F = -K*sin(θ)θ̂/|r|
    force = -K/d*(v - rdotv*rhat)

    return torque, force
end


"""
    torsion_torque(v, rhat, d, K)

Calculate the torque due to the bond angle torsion between two particles i and j, given the bond direction and the two orientations. Force will always be 0
"""
@inline @fastmath function torsion_torque(rhat, d, zi, zj, K_torque)
    rdotv = dot(zi,zj)
    rdotv = clamp(rdotv, -1, 1)

    

    # τ = K*sin(ϕ)*n̂
    torque = K_torque*(sqrt(1.0-rdotv^2))*rhat


    #print(norm(torque), "\n")

    #prevent numerical instability
    if (norm(torque)/K_torque < 1e-5)
        return zeros(3)
    end

    return torque
end


function norm_and_mag(r)
    d = norm(r)
    rhat = r ./ d
    return rhat, d
end

"""
    spring_force(rhat, d, l0::Real, k::Real)

Return force due to displacement r of a spring

# Arguments
- `rhat::BeadPos`: 3D direction vector of bond
- `d::Real`: length of bond
- `l0::Real`: rest length of spring
- `k::Real`: spring stiffness

# Returns
- `MVector{3, Float64}`: directed force
"""
@fastmath spring_force(rhat, d::Real, l0::Real, k::Real) = k*(d - l0) * rhat
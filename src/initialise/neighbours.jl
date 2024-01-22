
"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the lateral direction for bead `idx` of `total`
"""
function lateral_nn(idx::Int, total::Int; N::Int=13, S::Int=3)::Tuple{Int,Int}
    if idx % N == 0 
        if idx > total-(S-1)*N-1
            return (0, idx-1)
        end
        return (idx+1+(S-1)*N, idx-1)
    elseif idx % N == 1
        if idx < (S-1)*N+1
            return (idx+1, 0)
        end
        return (idx+1, idx-1-(S-1)*N)
    end
    return (idx+1, idx-1)
end

"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the longitudinal direction for bead `idx` of `total`
"""
function long_nn(idx::Int, total::Int; N::Int=13)::Tuple{Int,Int}
    if idx<N+1
        return idx+N, 0
    elseif idx>total-N
        return 0, idx-N
    end
    return idx+N, idx-N
end

"""
    lateral_nn(idx::Int, total::Int)::Tuple{Int,Int}

Return neighbours in the lateral direction for bead `idx` of `total` in a patch
"""
function lateral_nn_patch(idx::Int, N_lat::Int)::Tuple{Int,Int}
    if idx % N_lat == 0 
        return (0, idx-1)
    elseif idx % N_lat == 1
        return (idx+1, 0)
    end
    return (idx+1, idx-1)
end


"""
    long_nn_patch(idx::Int, N_lat::Int, N_long::Int)::Tuple{Int,Int}

Return neighbours in the long direction for bead `idx` of `total` in a patch
"""
function long_nn_patch(idx::Int, N_lat::Int, N_long::Int)::Tuple{Int,Int}
    if idx<N_lat+1
        return idx+N_lat, 0
    elseif idx>N_lat*(N_long-1)
        return 0, idx-N_lat
    end
    return idx+N_lat, idx-N_lat
end


"""
    neighbours(idx::Int, total::Int)

Return neighbours in the lateral and longitudinal direction for bead `idx` of `total`.
"""
function neighbours(idx::Int, total::Int, N::Int=13, S::Int=3)::Tuple{Tuple{Int,Int}, Tuple{Int, Int}}
    return lateral_nn(idx, total, N=N, S=S), long_nn(idx, total, N=N)
end


function intra_dimer_index(idx, total, α)
    lat, long = MicrotubuleSpringModel.neighbours(idx, total)
    return α ? 2+(lat[2]!=0)-(long[1]==0) : 1
end


function bond_indices(idx, total, α)
    lat, long = MicrotubuleSpringModel.neighbours(idx, total)
    bonds = []
    if long[1] != 0
        push!(bonds, α ? 3 : 2)
    end
    if lat[1] != 0
        push!(bonds, 1)
    end
    if long[2] != 0
        push!(bonds, α ? 2 : 3)
    end
    if lat[2] != 0
        push!(bonds, 1)
    end
    return return bonds
end



function iterateSDE!(lattice, consts, noise, dW, dt; S=3, N=13)
    rand!(noise, dW)
    @threads for i in 1:lastindex(lattice)
        if i > N
            F = linear_forces(lattice[i], lattice, consts)
            lattice[i].x .+= - F .* dt + dW[:,i]
        end
    end
end
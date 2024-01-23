#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

k_ab = 3.5
k_bb = 1.0/0.4
k_a = 3.6
k_b = 3.8
K_ab = 1/(1/(2*k_a)+1/(2*k_b)+1/k_ab)

K_bb = 1/(1/k_b+1/k_bb)

##############################################

# scales are all in nN, nm, ns 
# so energies will be nN.nm
# need to convert from that to Joules or k_B T etc

R = 2.0
η = 0.2
γ_r_inv = 1/(6π*η*R)
γ_θ_inv = 1/(8π*η*R^3)

@unpack S, N, dx, a = conf.lattice
r = S*a/N
R = N*dx/2π
δx = 2*R*sin(dx/(2*R))
l0_lat = sqrt(r^2 + δx^2)
ϕ = atan(r,δx)


Tk_B = 1.381e-5*300
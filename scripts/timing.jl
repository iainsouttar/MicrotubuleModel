#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end



#torsion bool is 0 if no torsion, 1 if is.
function record_times(torsion_bool, N, Nt)
    if torsion_bool == true
        conf = from_toml(MicrotubuleConfig, "config/eulerMT.toml")
        print("torsion")
    else
        conf = from_toml(MicrotubuleConfig, "config/eulerMTnoTorsion.toml")
        print("no Torsion")
    end
    print(conf.spring_consts)
    time = @elapsed begin
        for i in 1:N
            lattice, bead_info = initialise(conf)
            for j in 1:Nt
                iterate!(lattice, bead_info, conf, conf.iter_pars)
            end
        end
    end

    return time/N
end


times = zeros(15)
times1 = zeros(15)
Nts = 10 .^ range(1,4, 15)
for k in 1:15
    Nt = Nts[k]
    print("\n",Nt, "\n")
    times[k] = record_times(true, 200, Nt)
    times1[k] = record_times(false,200, Nt)

end

CairoMakie.activate!()
f = Figure(resolution=(1000,600))
ax = Axis(f[1,1])
lines!(ax, 10 .^ range(1,4, 15), (times), label = "With torsion")
lines!(ax, 10 .^ range(1,4, 15), (times1), label = "Without torsion")
ax.ylabel = "Time taken"
ax.xlabel = "Number of iterations"
axislegend(ax, position=:rb)

f
save("plots/timings.png", f)

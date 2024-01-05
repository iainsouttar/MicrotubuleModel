

function microtubule_length(beads, consts)
    N = consts.N
    Ntot = lastindex(beads)
    last = @view beads[Ntot-N:end]
    #tot = mapreduce((a,b)->norm(a.x-b.x), +, beads[1:N], last)
    # average axial distance between first and last rings of beads
    tot = mapreduce((a,b)->abs(a.x[3]-b.x[3]), +, beads[1:N], last)
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

function save_to_csv(filename::String, df::DataFrame; path="results/processed", header=true)
    CSV.write(path*"/"*filename, df, header=header)
end

function save_to_csv(filename::String, df::DataFrame, x::Vector; path="results/processed", header=true)
    df[!, "x"] = collect(x)
    select!(df, :x, Not([:x]))
    CSV.write(path*"/"*filename, df, header=header)
end


function save_to_csv(filename::String, x::Vector, ys; path="results/processed")
    data = hcat(x, ys...)
    writedlm(path*"/"*filename, data, ',')
end

function load_from_csv(filename::String; path="results/processed")
    return readdlm(path*"/"*filename)
end

save_params(fname, config; pth="results/processed") = to_toml(pth*"/"*fname, config)
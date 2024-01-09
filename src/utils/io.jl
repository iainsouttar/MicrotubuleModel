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
function save_to_csv(filename::String, df::DataFrame; path="results/processed", header=true, append=false)
    CSV.write(path*"/"*filename, df, header=header, append=append)
end

function save_to_csv(filename::String, df::DataFrame, x::Vector; path="results/processed", header=true)
    df[!, "x"] = collect(x)
    select!(df, :x, Not([:x]))
    CSV.write(path*"/"*filename, df, header=header)
end

function append_to_csv(filename::String, x::AbstractVector{T}, path="results/raw") where T
    data = Matrix{T}(zeros(T, (1,length(x))))
    for i in 1:lastindex(x)
        data[1,i] = x[i]
    end
    open(path*"/"*filename, "a") do io
        writedlm(io, data, ',')
    end
end

function append_to_csv(filename::String, x::AbstractVector{BeadPos}, path="results/raw")
    data = Matrix{Float64}(zeros(Float64, (1,length(x)*3)))
    for i in 1:lastindex(x)
        data[1,3*(i-1)+1:3*i] .= x[i]
    end
    open(path*"/"*filename, "a") do io
        writedlm(io, data, ',')
    end
end

function save_to_csv(filename::String, x::Vector, ys; path="results/processed")
    data = hcat(x, ys...)
    writedlm(path*"/"*filename, data, ',')
end

function load_from_csv(filename::String; path="results/processed")
    return readdlm(path*"/"*filename)
end

save_params(fname, config; pth="results/processed") = to_toml(pth*"/"*fname, config)
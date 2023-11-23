#  VS Code workaround: https://github.com/julia-vscode/julia-vscode/issues/800
if isdefined(@__MODULE__, :LanguageServer)
    @info "Using VS Code workaround..."
    include("../src/MicrotubuleSpringModel.jl")
    using .MicrotubuleSpringModel
else
    using MicrotubuleSpringModel
end

using Test
using LinearAlgebra

const tests = [
    "lattice",
    "orientation"
]


@testset "MicrotubuleSpringModel.jl" begin
    @testset "Test $t" for t in tests
        include("$t.jl")
    end
end
